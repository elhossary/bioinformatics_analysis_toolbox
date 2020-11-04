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
from wiggletools.wiggle_matrix import WiggleMatrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--merge_range", default=20, help="", type=int)
    parser.add_argument("--min_len", default=35, help="", type=int)
    parser.add_argument("--max_len", default=350, help="", type=int)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    gff_df["start"] = gff_df["start"].astype(int)
    chrom_sizes = get_chrom_sizes([os.path.abspath(args.refseq_in)])
    seqid_list = [x["seqid"] for x in chrom_sizes]
    wiggle_pathes = []
    for item in args.wigs_in:
        for sub_item in glob.glob(item):
            wiggle_pathes.append(os.path.abspath(sub_item))
    wiggles_parsed = [Wiggle(wiggle_path, chrom_sizes).get_wiggle(is_len_extended=True)
                      for wiggle_path in wiggle_pathes]
    wiggle_matrix = WiggleMatrix(wiggles_parsed, chrom_sizes, 4).wiggle_matrix_df
    f_scores_columns = [i for i in wiggle_matrix.columns.tolist() if "_forward" in i]
    r_scores_columns = [i for i in wiggle_matrix.columns.tolist() if "_reverse" in i]
    seqid_list = [i for i in seqid_list if i in wiggle_matrix["seqid"].unique().tolist()]
    seqid_list = [i for i in seqid_list if i in wiggle_matrix["seqid"].unique().tolist()]
    wiggle_matrix[r_scores_columns] = wiggle_matrix[r_scores_columns].abs()
    gff_df_len = gff_df.shape[0]
    slices_gff_df = pd.DataFrame(columns=col_names)
    anno_counter = 0
    for seqid in seqid_list:
        f_wig_df_slice = wiggle_matrix[wiggle_matrix["seqid"] == seqid].loc[:, f_scores_columns + ["location"]]
        r_wig_df_slice = wiggle_matrix[wiggle_matrix["seqid"] == seqid].loc[:, r_scores_columns + ["location"]]

        for idx in gff_df[gff_df["seqid"] == seqid].index:
            sys.stdout.flush()
            sys.stdout.write("\r" + f"Sequence ID {seqid} progress: {round(idx / gff_df_len * 100, 1)}%")
            start = gff_df.at[idx, "start"]
            end = gff_df.at[idx, "end"]
            strand = gff_df.at[idx, "strand"]
            list_out = []
            try:
                if strand == "+":
                    for score_column in f_scores_columns:
                        ret_result = slice_annotation_recursively(
                            f_wig_df_slice[f_wig_df_slice["location"]
                                .between(start, end)].loc[:, ["location", score_column]],
                            score_column, args.min_len, args.max_len)
                        if ret_result is not None and ret_result:
                            list_out.extend(ret_result)
                elif strand == "-":
                    for score_column in r_scores_columns:
                        ret_result = slice_annotation_recursively(
                            r_wig_df_slice[r_wig_df_slice["location"]
                                .between(start, end)].loc[:, ["location", score_column]],
                            score_column, args.min_len, args.max_len)
                        if ret_result is not None and ret_result:
                            list_out.extend(ret_result)
                else:
                    continue
            except Exception as e:
                print(e)
                continue
            slices_counter = 0
            if list_out is not None:
                list_out = sorted(list_out)
                list_out, _ = _merge_interval_lists(list_out, args.merge_range)
                list_out.extend(_)
                for i in list_out:
                    slices_counter += 1
                    anno_counter += 1
                    original_atrr = parse_attributes(gff_df.at[idx, 'attributes'])
                    attr = f"ID={seqid}_{strand}_{gff_df.at[idx, 'type']}_{anno_counter}" \
                           f";Name={original_atrr['name']}_slice_{slices_counter}" \
                           f";seq_len={i[1] - i[0] + 1}"
                    slices_gff_df = slices_gff_df.append(
                        {"seqid": seqid, "source": "ANNO_SLAYER", "type": gff_df.at[idx, 'type'],
                                   "start": i[0], "end": i[1], "score": ".", "strand": strand, "phase": ".",
                                   "attributes": attr}, ignore_index=True)
    slices_gff_df.to_csv(os.path.abspath(f"{args.gff_out}"), sep="\t", header=False, index=False)


def slice_annotation_recursively(coverage_df, score_col, min_len, max_len, ret_pos=None):
    if ret_pos is None:
        ret_pos = []
    max_score = coverage_df[score_col].max()
    cov_lst = coverage_df[score_col].tolist()
    peaks, peaks_prop = find_peaks(cov_lst, height=(max_score, max_score), width=(min_len, max_len), prominence=0)
    if len(cov_lst) < min_len or peaks.size == 0:
        return ret_pos
    width_heights = peaks_prop['width_heights'][0]
    peak_prominence = peaks_prop['prominences'][0]
    remove_limit = peaks_prop['width_heights'][0]
    max_score_loc = coverage_df.iloc[peaks[0]]['location']
    tmp_df = coverage_df[coverage_df[score_col] >= remove_limit]
    if tmp_df.shape[0] == coverage_df.shape[0]:
        return ret_pos
    single_sites = []
    for cg in consecutive_groups(tmp_df['location'].tolist()):
        consecutive_loc = list(cg)
        if len(consecutive_loc) == 1:
            single_sites.extend(consecutive_loc)
            continue
        if max_score_loc in consecutive_loc:
            cg_coverage = tmp_df[tmp_df["location"].isin(consecutive_loc)][score_col].tolist()
            mean_wid_prom = mean([peak_prominence, width_heights])
            mean_cg_ave_cov_prom = mean([peak_prominence, mean(cg_coverage)])
            #len_factor = peaks_prop['widths'][0] / 100
            tmp_df = tmp_df[tmp_df["location"].isin(consecutive_loc)]
            tmp_df = tmp_df[tmp_df[score_col] >= min([mean_wid_prom, mean_cg_ave_cov_prom])]
            new_cg = consecutive_groups(tmp_df['location'].tolist())
            new_cg = [list(cg) for cg in new_cg]
            for i in new_cg:
                if min_len <= max(i) - min(i) + 1 <= max_len:
                    ret_pos.append([min(i), max(i)])
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
        slice_annotation_recursively(coverage_df, score_col, min_len, max_len, ret_pos)
    return ret_pos


def _merge_interval_lists(list_in, merge_range):
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


def build_wiggle_matrix(pathes, chrom_sizes):

    pool = mp.Pool(processes=4)
    processes = []
    for wiggle_path in pathes:
        processes.append(pool.apply_async(
            Wiggle(os.path.abspath(wiggle_path), chrom_sizes).get_wiggle(),
            (os.path.abspath(wiggle_path), chrom_sizes, )))
    wigs_parsed = [p.get() for p in processes]
    mat = WiggleMatrix(wigs_parsed, chrom_sizes, 4).build_matrix()
    return mat.f_wiggle_matrix_df, mat.r_wiggle_matrix_df

def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()
exit(0)