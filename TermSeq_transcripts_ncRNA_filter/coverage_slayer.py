import pandas as pd
import os
from Bio import SeqIO
import sys
from more_itertools import consecutive_groups
from scipy.signal import find_peaks
import glob
from wiggletools.wiggle import Wiggle
import numpy as np
from functools import reduce
import argparse
import multiprocessing as mp


def main():
    np.set_printoptions(suppress=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--merge_range", default=0, help="", type=int)
    parser.add_argument("--min_len", default=50, help="", type=int)
    parser.add_argument("--max_len", default=350, help="", type=int)
    parser.add_argument("--threads", default=1, help="", type=int)
    parser.add_argument("--annotation_type", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    chrom_sizes = get_chrom_sizes([os.path.abspath(args.refseq_in)])
    seqid_list = [x["seqid"] for x in chrom_sizes]
    print("Parsing wiggles")
    wig_pool = mp.Pool(processes=args.threads)
    processes = []
    for item in args.wigs_in:
        for sub_item in glob.glob(item):
            processes.append(wig_pool.apply_async(create_wiggle_obj,
                                                  args=(os.path.abspath(sub_item), chrom_sizes)))
    wiggles_parsed = [p.get() for p in processes]
    wig_pool.close()
    print("Differentiating wiggles")
    f_wiggles = [i for i in wiggles_parsed if i.orientation == "f"]
    r_wiggles = [i for i in wiggles_parsed if i.orientation == "r"]
    del wiggles_parsed
    print("Generating forward arrays")
    f_arrays = [convert_wiggle_obj_to_arr(i, args) for i in f_wiggles]
    print("Generating reverse arrays")
    r_arrays = [convert_wiggle_obj_to_arr(i, args) for i in r_wiggles]
    f_raw_pos = {}
    r_raw_pos = {}
    pool = mp.Pool(processes=args.threads)
    f_proc = []
    r_proc = []
    for seqid in seqid_list:
        if not seqid in f_raw_pos.keys():
            f_raw_pos[seqid] = []
        print(f"Predicting peaks for {seqid}+")
        for arr in f_arrays:
            for loc in arr[1][seqid]:
                f_proc.append(pool.apply_async(predict_locs_arr_slice,
                                               args=(arr[0][seqid], loc[0], loc[1], args)))
            proc_len = len(f_proc)
            counter = 0
            for p in f_proc:
                counter += 1
                sys.stdout.flush()
                sys.stdout.write("\r" + f"Progress: {round(counter / proc_len * 100, 1)}%")
                f_raw_pos[seqid].extend(p.get())

    for seqid in seqid_list:
        if not seqid in r_raw_pos.keys():
            r_raw_pos[seqid] = []
        print(f"Predicting peaks for {seqid}-")
        for arr in r_arrays:
            for loc in arr[1][seqid]:
                r_proc.append(pool.apply_async(predict_locs_arr_slice,
                                               args=(arr[0][seqid], loc[0], loc[1], args)))
            proc_len = len(r_proc)
            counter = 0
            for p in r_proc:
                counter += 1
                sys.stdout.flush()
                sys.stdout.write("\r" + f"Progress: {round(counter / proc_len * 100, 1)}%")
                r_raw_pos[seqid].extend(p.get())

    pool.close()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    peaks_gff_df = pd.DataFrame(columns=col_names)
    print("Generating annotations from merged overlapping positions for forward")
    for k, v in f_raw_pos.items():
        peaks_gff_df = peaks_gff_df.append(
            generate_annotations_from_positions(v, k, "+", args.annotation_type, 0),
            ignore_index=True)
    print("Generating annotations from merged overlapping positions for reverse")
    for k, v in r_raw_pos.items():
        peaks_gff_df = peaks_gff_df.append(
            generate_annotations_from_positions(v, k, "-", args.annotation_type, 0),
            ignore_index=True)
    print("Exporting to GFF")
    peaks_gff_df["len"] = peaks_gff_df["end"] - peaks_gff_df["start"] + 1
    len_range = range(args.min_len, args.max_len + 1, 1)
    peaks_gff_df = peaks_gff_df[peaks_gff_df["len"].isin(len_range)]
    peaks_gff_df.drop("len", axis=1, inplace=True)
    peaks_gff_df.sort_values(["seqid", "start", "end"], inplace=True)
    peaks_gff_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)


def predict_locs_arr_slice(arr, start, end, args):
    arr_slice = arr[np.logical_and(arr[:, 0] >= start, arr[:, 0] <= end)].copy()
    return generate_locs(arr_slice, args)


def generate_locs(wig_arr_slice, args):
    dist = args.min_len * 2
    length_range = range(args.min_len, args.max_len + 1, 1)
    ret_locs = []
    # Assuming that the incoming array have the raw coverage values at column 1,
    # where the rising step height at 2, and falling step height at 3
    # of course column 0 is for location
    if len(set(wig_arr_slice[:, 1].tolist())) > 1:
        wig_peaks, wig_peaks_prop = find_peaks(wig_arr_slice[:, 1], prominence=0, distance=dist)
        rising_peaks, rising_peaks_prop = find_peaks(wig_arr_slice[:, 2], prominence=0, distance=dist)
        falling_peaks, falling_peaks_prop = find_peaks(wig_arr_slice[:, 3], prominence=0, distance=dist)
        wig_locs = [wig_arr_slice[i, 0] for i in wig_peaks]
        rising_locs = [wig_arr_slice[i, 0] for i in rising_peaks]
        falling_locs = [wig_arr_slice[i, 0] for i in falling_peaks]
        for l in wig_locs:
            uppers = [i for i in falling_locs if l <= i]
            lowers = [i for i in rising_locs if l >= i]
            if not uppers or not lowers:
                continue
            lower_loc = max(lowers)
            upper_loc = min(uppers)
            if upper_loc - lower_loc + 1 in length_range:
                ret_locs.append([lower_loc, upper_loc])
    elif len(set(wig_arr_slice[:, 1].tolist())) == 1:
        locs = wig_arr_slice[:, 0].tolist()
        lower_loc = min(locs)
        upper_loc = max(locs)
        if upper_loc - lower_loc + 1 in length_range:
            ret_locs.append([min(locs), max(locs)])
    else:
        print("Fatal Error")
    return ret_locs

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
                 "source": "COVERAGE_SLAYER",
                 "type": new_type,
                 "start": i[0],
                 "end": i[1],
                 "score": ".",
                 "strand": strand,
                 "phase": ".",
                 "attributes": attr})
    return anno_list


def convert_wiggle_obj_to_arr(wig_obj, args):
    arr_dict = {}
    wig_cols = ["variableStep_chrom", "location", "score"]
    wig_df = wig_obj.get_wiggle().loc[:, wig_cols]
    wig_df["score"] = wig_df["score"].abs()
    merged_df = reduce(lambda x, y: pd.merge(x, y, on=["variableStep_chrom", "location"], how='left'),
                       [wig_df.loc[:, wig_cols],
                        wig_obj.to_step_height(3, "start_end").loc[:, wig_cols],
                        wig_obj.to_step_height(3, "end_start").loc[:, wig_cols]])
    for seqid in merged_df["variableStep_chrom"].unique():
        tmp = merged_df[merged_df["variableStep_chrom"] == seqid].drop("variableStep_chrom", axis=1)
        arr_dict[seqid] = tmp.to_numpy(copy=True)
    return arr_dict, generate_locations(wig_df, args)


def generate_locations(wig_df, args):
    stretches_locs = {}
    wig_df_minified = wig_df[wig_df["score"] > 0.0].loc[:, ["variableStep_chrom", "location", "score"]]
    for seqid in wig_df_minified["variableStep_chrom"].unique().tolist():
        cons_grps = consecutive_groups(
            wig_df_minified[wig_df_minified["variableStep_chrom"] == seqid]["location"].tolist())
        stretches_locs[seqid] = []
        for cg in cons_grps:
            cg_l = list(cg)
            if len(cg_l) >= args.min_len:
                stretches_locs[seqid].append([min(cg_l), max(cg_l)])
    return stretches_locs



def create_wiggle_obj(fpath, chrom_sizes):
    return Wiggle(fpath, chrom_sizes, is_len_extended=True)


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


main()