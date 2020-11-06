import argparse
import pandas as pd
import os
import multiprocessing as mp
from Bio import SeqIO
import sys
from statistics import mean, median
from more_itertools import consecutive_groups
from scipy.signal import find_peaks
import glob
from wiggletools.wiggle import Wiggle
from itertools import product



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--merge_range", default=20, help="", type=int)
    parser.add_argument("--min_len", default=50, help="", type=int)
    parser.add_argument("--max_len", default=350, help="", type=int)
    parser.add_argument("--threads", default=1, help="", type=int)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    gff_df["start"] = gff_df["start"].astype(int)
    gff_df["end"] = gff_df["end"].astype(int)
    chrom_sizes = get_chrom_sizes([os.path.abspath(args.refseq_in)])
    seqid_list = gff_df["seqid"].unique().tolist()
    wig_pool = mp.Pool(processes=args.threads)
    processes = []
    for item in args.wigs_in:
        for sub_item in glob.glob(item):
            processes.append(wig_pool.apply_async(create_wiggle_obj, args=(os.path.abspath(sub_item), chrom_sizes)))
    wiggles_parsed = [p.get() for p in processes]
    wiggles_parsed = [df[df["variableStep_chrom"].isin(seqid_list)].copy() for df in wiggles_parsed]
    f_wiggles = [x for x in wiggles_parsed if any(x["score"] > 0)]
    r_wiggles = [x for x in wiggles_parsed if any(x["score"] < 0)]
    for wig in r_wiggles:
        wig["score"] = wig["score"].abs()
    slices_gff_df = pd.DataFrame(columns=col_names)
    anno_counter = 0
    combinations = product(seqid_list, ["+", "-"])
    for comb in combinations:
        print(f"Trimming for sequence ID {comb[0]}{comb[1]}")
        wigs_selection = f_wiggles if comb[1] == "+" else r_wiggles
        wiggles_slices = [df[df["variableStep_chrom"] == comb[0]] for df in wigs_selection]
        anno_slice = gff_df[(gff_df["seqid"] == comb[0]) & (gff_df["strand"] == comb[1])]
        gff_len = max(anno_slice.index.tolist())
        slicer_pool = mp.Pool(processes=args.threads)
        slicer_processes = []
        for idx in anno_slice.index:
            sys.stdout.flush()
            sys.stdout.write("\r" + f"Progress: {round(idx / gff_len * 100, 1)}%")
            anno_counter += 1
            row = gff_df.loc[idx, :]
            slicer_processes.append(slicer_pool.apply_async(slicer, (row, wiggles_slices, args, comb, anno_counter)))
            slicer_res = [p.wait() for p in slicer_processes]
            slices_gff_df.append(slicer_res, ignore_index=True)

    slices_gff_df.to_csv(os.path.abspath(f"{args.gff_out}"), sep="\t", header=False, index=False)


def slicer(row, wiggles_slices, args, comb, anno_counter):
    start = row["start"]
    end = row["end"]
    strand = row["strand"]
    list_out = []
    for wig in wiggles_slices:
        wig_pos_slice = wig[wig["location"].between(start, end)].loc[:, ["location", "score"]]
        if wig_pos_slice.empty:
            continue
        try:
            ret_result = slice_annotation_recursively(wig_pos_slice, args.min_len, args.max_len)
            if ret_result is not None and ret_result:
                list_out.extend(ret_result)
            else:
                continue
        except Exception as e:
            print(f"Warning: {e}")
            continue
    slices_counter = 0
    if list_out:
        list_out, _ = _merge_interval_lists(list_out, args.merge_range)
        list_out.extend(_)
        for i in list_out:
            slices_counter += 1
            original_atrr = parse_attributes(row['attributes'])
            attr = f"ID={comb[1]}_{strand}_{row['type']}_{anno_counter}" \
                   f";Name={original_atrr['name']}_slice_{slices_counter}" \
                   f";seq_len={i[1] - i[0] + 1}"
            return \
                {"seqid": comb[1],
                 "source": "ANNO_SLAYER",
                 "type": row['type'],
                 "score": ".",
                 "strand": strand,
                 "phase": ".",
                 "attributes": attr}


def create_wiggle_obj(fpath, chrom_sizes):
    return Wiggle(fpath, chrom_sizes).get_wiggle(is_len_extended=True)


def slice_annotation_recursively(coverage_df, min_len, max_len, ret_pos=None):
    if ret_pos is None:
        ret_pos = []
    max_score = coverage_df["score"].max()
    cov_lst = coverage_df["score"].tolist()
    peaks, peaks_prop = find_peaks(cov_lst, height=(max_score, max_score), width=(min_len, max_len), prominence=0)
    if len(cov_lst) < min_len or peaks.size == 0:
        return ret_pos
    width_heights = peaks_prop['width_heights'][0]
    peak_prominence = peaks_prop['prominences'][0]
    max_score_loc = int(coverage_df.iloc[peaks[0]]['location'])
    tmp_df = coverage_df[coverage_df["score"] >= width_heights]
    if tmp_df.shape[0] == coverage_df.shape[0]:
        return ret_pos
    single_sites = []
    for cg in consecutive_groups(tmp_df['location'].tolist()):
        consecutive_loc = list(cg)
        if len(consecutive_loc) == 1:
            single_sites.extend(consecutive_loc)
            continue
        if max_score_loc in consecutive_loc:
            cg_coverage = tmp_df[tmp_df["location"].isin(consecutive_loc)]["score"].tolist()
            mean_wid_prom = mean([peak_prominence, width_heights])
            mean_cg_ave_cov_prom = mean([peak_prominence, mean(cg_coverage)])
            len_factor = peaks_prop['widths'][0] / 100
            tmp_df = tmp_df[tmp_df["location"].isin(consecutive_loc)]
            tmp_df = tmp_df[tmp_df["score"] >= min([mean_wid_prom, mean_cg_ave_cov_prom]) * len_factor]
            new_cg = consecutive_groups(tmp_df['location'].tolist())
            new_cg = [list(ncg) for ncg in new_cg]
            for ncg in new_cg:
                if min_len <= max(ncg) - min(ncg) + 1 <= max_len:
                    ret_pos.append([min(ncg), max(ncg)])
            """
            print(new_cg)
            new_cg = [i for i in new_cg if min_len <= max(i) - min(i) + 1 <= max_len]
            print(len(new_cg))
            x = [tmp_df["location"].min(), tmp_df["location"].max()]
            if min_len <= x[-1] - x[0] + 1 <= max_len:
                ret_pos.append(x)
            else:
                if x[-1] - x[0] + 1 > max_len:
                    print("warning too long")
            """
            single_sites.extend(consecutive_loc)
            break
    coverage_df = coverage_df[~coverage_df["location"].isin(single_sites)].copy()
    if not coverage_df.empty:
        slice_annotation_recursively(coverage_df, min_len, max_len, ret_pos)
    return ret_pos


def _merge_interval_lists(list_in, merge_range):
    list_in = sorted(list_in)
    merge_range += 2
    list_out = []
    overlap_indices = []
    for loc in list_in:
        if len(list_out) == 0:
            list_out.append(loc)
        else:
            if loc[0] in range(list_out[-1][0], list_out[-1][-1] + merge_range):
                list_out[-1][-1] = max([loc[-1], list_out[-1][-1]])
                overlap_indices.append(list_out.index(list_out[-1]))
            else:
                list_out.append(loc)
    overlap_indices = list(set(overlap_indices))
    overlap_indices.sort()
    overlaps_list_out = [list_out[i] for i in overlap_indices]
    return overlaps_list_out, [i for i in list_out if i not in overlaps_list_out]

def get_chrom_sizes(fasta_pathes):
    ret_list = []
    for fasta_path in fasta_pathes:
        print(f"==> Parsing reference sequence: {os.path.basename(fasta_path)}")
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            ret_list.append({"seqid": seq_record.id,
                             "size": len(seq_record.seq),
                             "fasta": os.path.basename(fasta_path)})
    return ret_list


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()
print("\n")
exit(0)