import pandas as pd
import os
from Bio import SeqIO
from scipy.signal import find_peaks
import glob
from wiggletools.wiggle import Wiggle
import numpy as np
from functools import reduce
import argparse
import multiprocessing as mp
from statistics import mean


def main():
    np.set_printoptions(suppress=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--merge_range", default=0, help="", type=int)
    parser.add_argument("--min_len", default=35, help="", type=int)
    parser.add_argument("--max_len", default=300, help="", type=int)
    parser.add_argument("--peak_distance", default=70, help="", type=int)
    parser.add_argument("--ignore_coverage", default=5, help="", type=int)
    parser.add_argument("--threads", default=1, help="", type=int)
    parser.add_argument("--annotation_type", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    chrom_sizes = get_chrom_sizes([os.path.abspath(args.refseq_in)])

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
    del f_wiggles
    print("Generating reverse arrays")
    r_arrays = [convert_wiggle_obj_to_arr(i, args) for i in r_wiggles]
    del r_wiggles
    f_raw_predicted_locs = {}
    for arr_obj in f_arrays:
        print(f"Predicting peaks for condition: '{arr_obj[1]}'")
        for seqid_key in arr_obj[0].keys():
            print(f"\tProcessing sequence ID: '{seqid_key}'")
            if seqid_key not in f_raw_predicted_locs.keys():
                f_raw_predicted_locs[seqid_key] = []
            f_raw_predicted_locs[seqid_key].extend(generate_locs(arr_obj[0][seqid_key], args, False, arr_obj[1]))
    r_raw_predicted_locs = {}
    for arr_obj in r_arrays:
        print(f"Predicting peaks for condition: '{arr_obj[1]}'")
        for seqid_key in arr_obj[0].keys():
            print(f"\tProcessing sequence ID: '{seqid_key}'")
            if seqid_key not in r_raw_predicted_locs.keys():
                r_raw_predicted_locs[seqid_key] = []
            r_raw_predicted_locs[seqid_key].extend(generate_locs(arr_obj[0][seqid_key], args, True, arr_obj[1]))

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    peaks_gff_df = pd.DataFrame(columns=col_names)
    for k, v in f_raw_predicted_locs.items():
        print(f"Generating annotations from merged overlapping positions for forward sequence ID: {k}")
        peaks_gff_df = peaks_gff_df.append(
            generate_annotations_from_positions(v, k, "+", args.annotation_type, args),
            ignore_index=True)
    for k, v in r_raw_predicted_locs.items():
        print(f"Generating annotations from merged overlapping positions for reverse sequence ID: {k}")
        peaks_gff_df = peaks_gff_df.append(
            generate_annotations_from_positions(v, k, "-", args.annotation_type, args),
            ignore_index=True)
    print("Filtering by length")
    peaks_gff_df["len"] = peaks_gff_df["end"] - peaks_gff_df["start"] + 1
    len_range = range(args.min_len, args.max_len + 1, 1)
    peaks_gff_df = peaks_gff_df[peaks_gff_df["len"].isin(len_range)]
    peaks_gff_df.drop("len", axis=1, inplace=True)
    print("Sorting annotated peaks")
    peaks_gff_df.sort_values(["seqid", "start", "end"], inplace=True)
    print(f"Total {peaks_gff_df.shape[0]} peaks predicted, exporting to GFF")
    peaks_gff_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)


def generate_locs(coverage_array, args, is_reversed, cond_name):
    length_range = range(args.min_len, args.max_len + 1, 1)
    ret_locs = []
    # Assuming that the incoming array have the raw coverage values at column 1,
    # where the rising step height at 2, and falling step height at 3
    # of course column 0 is for location
    wig_peaks, wig_peaks_prop =\
        find_peaks(coverage_array[:, 1],
                   distance=args.peak_distance, width=(args.min_len, args.max_len), rel_height=1)
    rising_peaks, rising_peaks_prop =\
        find_peaks(coverage_array[:, 2],
                   distance=args.peak_distance, width=(1, 7), rel_height=1)
    falling_peaks, falling_peaks_prop =\
        find_peaks(coverage_array[:, 3],
                   distance=args.peak_distance, width=(1, 7), rel_height=1)
    wig_locs = [coverage_array[i, 0] for i in wig_peaks]
    rising_locs = [coverage_array[i, 0] for i in rising_peaks]
    falling_locs = [coverage_array[i, 0] for i in falling_peaks]
    for loc in wig_locs:
        lowers = [i for i in falling_locs if loc >= i] if is_reversed else [i for i in rising_locs if loc >= i]
        uppers = [i for i in rising_locs if loc <= i] if is_reversed else [i for i in falling_locs if loc <= i]
        if not uppers or not lowers:
            continue
        lower_loc = int(max(lowers))
        upper_loc = int(min(uppers))
        average_coverage = round(mean(
            coverage_array[np.logical_and(lower_loc <= coverage_array[:, 0], coverage_array[:, 0] <= upper_loc)]
            [:, 1].copy().tolist()), 2)
        if upper_loc - lower_loc + 1 in length_range:
            ret_locs.append([lower_loc, upper_loc, average_coverage,
                             f"{cond_name}_with_best_ave_coverage_{average_coverage}"])
    return ret_locs


def generate_annotations_from_positions(list_out, seqid, strand, new_type, args):
    strand_letter_func = lambda x: "F" if x == "+" else "R"
    anno_counter = 0
    anno_list = []
    if list_out:
        list_out, _ = _merge_interval_lists(list_out, args)
        list_out.extend(_)
        for i in list_out:
            anno_counter += 1
            attr = f"ID={seqid}{strand_letter_func(strand)}_{anno_counter}" \
                   f";Name={seqid}_{strand_letter_func(strand)}_{new_type}_{anno_counter}" \
                   f";seq_len={i[1] - i[0] + 1}" \
                   f";conditions={i[3]}"
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
    step_size = int(args.min_len / 2)
    arr_dict = {}
    wig_cols = ["variableStep_chrom", "location", "score"]
    wig_df = wig_obj.get_wiggle()
    cond_name = wig_df.iat[0, 1]
    wig_df = wig_df.loc[:, wig_cols]
    wig_df["score"] = wig_df["score"].abs()
    wig_df.loc[wig_df['score'] <= args.ignore_coverage, ["score"]] = 0.0
    merged_df = reduce(lambda x, y: pd.merge(x, y, on=["variableStep_chrom", "location"], how='left'),
                       [wig_df.loc[:, wig_cols],
                        wig_obj.to_step_height(step_size, "start_end").loc[:, wig_cols],
                        wig_obj.to_step_height(step_size, "end_start").loc[:, wig_cols]])
    merged_df["location"] = merged_df["location"].astype(int)
    for seqid in merged_df["variableStep_chrom"].unique():
        tmp = merged_df[merged_df["variableStep_chrom"] == seqid].drop("variableStep_chrom", axis=1)
        arr_dict[seqid] = np.absolute(tmp.to_numpy(copy=True))
    return arr_dict, cond_name


def create_wiggle_obj(fpath, chrom_sizes):
    return Wiggle(fpath, chrom_sizes, is_len_extended=True)


def _merge_interval_lists(list_in, args):
    # Indices of location information
    start = 0
    end = 1
    cov_mean = 2
    cond_name = 3
    last_loc = -1
    ################
    select_func = lambda x, y: x if x[cov_mean] > y[cov_mean] else y if y[cov_mean] > x[cov_mean] else None
    print("Merging overlapping locations")
    merge_range = args.merge_range + 1
    list_out = []
    list_in = sorted(list_in)
    overlap_indices = []
    for new_loc in list_in:
        if list_out:
            if new_loc[start] in range(list_out[last_loc][start], list_out[last_loc][end] + merge_range):
                selected_loc = select_func(list_out[last_loc], new_loc)
                if selected_loc is None:
                    list_out[last_loc][end] = max([new_loc[end], list_out[last_loc][end]])
                    list_out[last_loc][cond_name] += f",{new_loc[cond_name]}"
                else:
                    list_out[last_loc][start], list_out[last_loc][end], list_out[last_loc][cond_name] = \
                        selected_loc[start], selected_loc[end], selected_loc[cond_name]

                if list_out[last_loc] not in overlap_indices:
                    overlap_indices.append(list_out.index(list_out[last_loc]))
                overlap_indices.append(list_out.index(new_loc))
            else:
                # check neighbors if no overlap
                if new_loc[end] - list_out[last_loc][start] + 1 in range(0, args.max_len + 1, 1) and \
                        new_loc[start] - list_out[last_loc][end] + 1 in range(0, args.min_len):
                    selected_loc = select_func(list_out[last_loc], new_loc)

                    if selected_loc is None:
                        list_out.append(new_loc)
                    else:
                        list_out[last_loc][start], list_out[last_loc][end], list_out[last_loc][cond_name] = \
                            selected_loc[start], selected_loc[end], selected_loc[cond_name]
                else:
                    list_out.append(new_loc)
        else:
            list_out.append(new_loc)

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