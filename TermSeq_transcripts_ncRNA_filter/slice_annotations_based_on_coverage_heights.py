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
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--merge_range", default=20, help="", type=int)
    parser.add_argument("--min_len", default=50, help="", type=int)
    parser.add_argument("--max_len", default=350, help="", type=int)
    parser.add_argument("--threads", default=1, help="", type=int)
    parser.add_argument("--rename_type", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    gff_df["start"] = gff_df["start"].astype(int)
    chrom_sizes = get_chrom_sizes([os.path.abspath(args.refseq_in)])
    seqid_list = [x["seqid"] for x in chrom_sizes]
    wig_pool = mp.Pool(processes=args.threads)
    processes = []
    for item in args.wigs_in:
        for sub_item in glob.glob(item):
            processes.append(wig_pool.apply_async(create_wiggle_obj, args=(os.path.abspath(sub_item), chrom_sizes)))
    wiggles_parsed = [p.get() for p in processes]
    wiggle_matrix = WiggleMatrix(wiggles_parsed, chrom_sizes, args.threads).wiggle_matrix_df
    wig_pool.close()
    f_scores_columns = [i for i in wiggle_matrix.columns.tolist() if "_forward" in i]
    r_scores_columns = [i for i in wiggle_matrix.columns.tolist() if "_reverse" in i]
    seqid_list = [i for i in seqid_list if i in wiggle_matrix["seqid"].unique().tolist()]
    seqid_list = [i for i in seqid_list if i in wiggle_matrix["seqid"].unique().tolist()]
    wiggle_matrix[r_scores_columns] = wiggle_matrix[r_scores_columns].abs()
    slices_gff_df = pd.DataFrame(columns=col_names)
    slicer_pool = mp.Pool(processes=args.threads)
    slicer_processes = []
    print("Spawning processes")
    for seqid in seqid_list:
        f_wig_df_slice = wiggle_matrix[wiggle_matrix["seqid"] == seqid].loc[:, f_scores_columns + ["location"]]
        r_wig_df_slice = wiggle_matrix[wiggle_matrix["seqid"] == seqid].loc[:, r_scores_columns + ["location"]]
        for idx in gff_df[gff_df["seqid"] == seqid].index:
            slicer_processes.append(
                slicer_pool.apply_async(row_process,
                                        (gff_df.loc[idx],
                                         f_wig_df_slice, r_wig_df_slice,
                                         f_scores_columns, r_scores_columns,
                                         args)))
    proc_len = len(slicer_processes)
    counter = 0
    for p in slicer_processes:
        counter += 1
        sys.stdout.flush()
        sys.stdout.write("\r" + f"Progress: {round(counter / proc_len * 100, 1)}%")
        slices_gff_df = slices_gff_df.append(p.get())
    slicer_pool.close()
    slices_gff_df.sort_values(["seqid", "start", "end"], inplace=True)
    slices_gff_df.to_csv(os.path.abspath(f"{args.gff_out}"), sep="\t", header=False, index=False)


def row_process(row, f_wig_df_slice, r_wig_df_slice, f_scores_columns, r_scores_columns, args):
    start = row["start"]
    end = row["end"]
    strand = row["strand"]
    seqid = row["seqid"]
    wig_selection = f_wig_df_slice[f_wig_df_slice["location"].isin(range(start, end + 1))].copy() if strand == "+" \
        else r_wig_df_slice[r_wig_df_slice["location"].isin(range(start, end + 1))].copy()
    col_selection = f_scores_columns if strand == "+" else r_scores_columns
    list_out = []
    results_collection = [
        slice_annotation_recursively(wig_selection.loc[:, ["location", x]], x, args.min_len, args.max_len) for x in
        col_selection]
    for x in results_collection:
        list_out.extend(x)
    return generate_annotations_from_positions(list_out, seqid, strand, args.rename_type, args.merge_range)


def generate_annotations_from_positions(list_out, seqid, strand, new_type, merge_range):
    strand_letter_func = lambda x: "F" if x == "+" else "R"
    anno_counter = 0
    anno_list = []
    if list_out:
        list_out, _ = _merge_interval_lists(list_out, merge_range)
        list_out.extend(_)
        for i in list_out:
            anno_counter += 1
            attr = f"ID={seqid}{strand_letter_func(strand)}_{anno_counter}" \
                   f";Name={seqid}_{strand_letter_func(strand)}_{new_type}_{anno_counter}" \
                   f";seq_len={i[1] - i[0] + 1}"
            anno_list.append(
                {"seqid": seqid,
                 "source": "ANNO_SLAYER",
                 "type": new_type,
                 "start": i[0],
                 "end": i[1],
                 "score": ".",
                 "strand": strand,
                 "phase": ".",
                 "attributes": attr})
    return anno_list


def create_wiggle_obj(fpath, chrom_sizes):
    return Wiggle(fpath, chrom_sizes).get_wiggle(is_len_extended=True)


def slice_annotation_recursively(coverage_df, score_col, min_len, max_len, ret_pos=None):
    try:
        if ret_pos is None:
            ret_pos = []
        max_score = coverage_df[score_col].max()
        cov_lst = coverage_df[score_col].tolist()
        peaks, peaks_prop = find_peaks(cov_lst, height=(max_score, max_score), width=(min_len, max_len), prominence=0)
        if len(cov_lst) < min_len or peaks.size == 0:
            return ret_pos
        width_heights = peaks_prop['width_heights'][0]
        peak_prominence = peaks_prop['prominences'][0]
        max_score_loc = int(coverage_df.iloc[peaks[0]]['location'])
        tmp_df = coverage_df[coverage_df[score_col] >= width_heights]
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
                len_factor = peaks_prop['widths'][0] / 100
                tmp_df = tmp_df[tmp_df["location"].isin(consecutive_loc)]
                tmp_df = tmp_df[tmp_df[score_col] >= min([mean_wid_prom, mean_cg_ave_cov_prom]) * len_factor]
                new_cg = consecutive_groups(tmp_df['location'].tolist())
                new_cg = [list(ncg) for ncg in new_cg]
                for ncg in new_cg:
                    if min_len <= max(ncg) - min(ncg) + 1 <= max_len:
                        ret_pos.append([min(ncg), max(ncg)])
                single_sites.extend(consecutive_loc)
                break
        coverage_df = coverage_df[~coverage_df["location"].isin(single_sites)].copy()
        if not coverage_df.empty:
            slice_annotation_recursively(coverage_df, score_col, min_len, max_len, ret_pos)
        return ret_pos
    except Exception as e:
        print(f"Warning: {e}")
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