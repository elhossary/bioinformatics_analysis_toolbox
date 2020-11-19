import pandas as pd
import os
from Bio import SeqIO
from scipy.signal import find_peaks
import glob
import sys
from wiggletools.wiggle import Wiggle
import numpy as np
from functools import reduce
import argparse
import multiprocessing as mp
from more_itertools import consecutive_groups
from statistics import mean


def main():
    np.set_printoptions(suppress=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--merge_range", default=0, help="", type=int)
    parser.add_argument("--max_gap", default=20, help="", type=int)
    parser.add_argument("--min_len", default=30, help="", type=int)
    parser.add_argument("--max_len", default=350, help="", type=int)
    parser.add_argument("--ignore_coverage", default=10, help="", type=int)
    parser.add_argument("--threads", default=1, help="", type=int)
    parser.add_argument("--annotation_type", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    args.merge_range += 1
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
            print("\n")
    r_raw_predicted_locs = {}
    for arr_obj in r_arrays:
        print(f"Predicting peaks for condition: '{arr_obj[1]}'")
        for seqid_key in arr_obj[0].keys():
            print(f"\tProcessing sequence ID: '{seqid_key}'")
            if seqid_key not in r_raw_predicted_locs.keys():
                r_raw_predicted_locs[seqid_key] = []
            r_raw_predicted_locs[seqid_key].extend(generate_locs(arr_obj[0][seqid_key], args, True, arr_obj[1]))
            print("\n")

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
    # Assuming that the incoming array have the raw coverage values at column 1,
    # where the rising step height at 2, and falling step height at 3
    # of course column 0 is for location
    mean_func = lambda pos: mean(coverage_array[pos[0] - 1:pos[1] - 1, 1].tolist())
    rising_peaks, _ = find_peaks(coverage_array[:, 2])
    falling_peaks, _ = find_peaks(coverage_array[:, 3])
    rising_peaks_list = rising_peaks.tolist()
    falling_peaks_list = falling_peaks.tolist()
    falling_peaks_set = set(falling_peaks_list)
    rising_peaks_len = len(rising_peaks_list)
    counter = 0
    possible_locations = []
    for rp in rising_peaks_list:
        counter += 1
        sys.stdout.flush()
        sys.stdout.write("\r" + f"\t\tProgress: {round(counter / rising_peaks_len * 100, 1)}%")
        rp_height = coverage_array[rp, 2]
        range_start = max(rp - args.max_len, 0) if is_reversed else rp + args.min_len
        range_end = rp - args.min_len if is_reversed else rp + args.max_len
        fp_range = set(range(range_start, range_end, 1))
        possible_fp = list(falling_peaks_set.intersection(fp_range))
        if not possible_fp:
            continue
        possible_fp_heights = coverage_array[possible_fp, 3].tolist()
        closest_falling_cov = min(possible_fp_heights, key=lambda x: abs(x - rp_height))
        closest_cov_fp = possible_fp[possible_fp_heights.index(closest_falling_cov)]
        upper_loc = int(coverage_array[rp, 0]) if is_reversed else int(coverage_array[closest_cov_fp, 0])
        lower_loc = int(coverage_array[closest_cov_fp, 0]) if is_reversed else int(coverage_array[rp, 0])
        possible_locations.append([lower_loc, upper_loc, mean_func([lower_loc, upper_loc]), cond_name])
    possible_locations = select_highest_mean_coverage_from_overlaps(possible_locations)
    possible_locations = merge_locations_with_gaps(possible_locations, args, coverage_array)
    print(f"\t\tLocations found: {len(possible_locations)}")
    return possible_locations


def select_highest_mean_coverage_from_overlaps(list_in):
    list_in = sorted(list_in, key=lambda l: (l[0], l[1]))
    start = 0
    end = 1
    cmean = 2
    last = -1
    select_func = lambda x, y: 1 if x[cmean] > y[cmean] else 2 if y[cmean] > x[cmean] else 0
    list_out = []
    for loc in list_in:
        if not list_out:
            # append first element
            list_out.append(loc)
            continue
        if loc[start] in range(list_out[last][start], list_out[last][end] + 2):
            selected = select_func(list_out[last], loc)
            if selected == 2:
                list_out[last] = loc
        else:
            list_out.append(loc)
    return list_out


def merge_locations_with_gaps(list_in, args, coverage_array):
    mean_func = lambda pos: mean(coverage_array[pos[0] - 1:pos[1] - 1, 1].tolist())
    gap_range = range(0, args.max_gap + 2, 1)
    length_range = range(args.min_len, args.max_len + 1, 1)
    list_in = sorted(list_in, key=lambda l: (l[0], l[1]))
    start = 0
    end = 1
    cmean = 2
    last = -1
    list_out = []
    for loc in list_in:
        if not list_out:
            # append first element
            list_out.append(loc)
            continue
        # check the gap between the new candidate and latest appended one
        if loc[start] - list_out[last][end] - 1 in gap_range and \
                loc[start] - list_out[last][0] + 1 in length_range:
            list_out[last][end] = loc[end]
            list_out[last][cmean] = mean_func([list_out[last][start], list_out[last][end]])
        else:
            list_out.append(loc)
    return list_out


def generate_annotations_from_positions(list_out, seqid, strand, new_type, args):
    strand_letter_func = lambda x: "F" if x == "+" else "R"
    anno_counter = 0
    anno_list = []
    if list_out:
        list_out, _ = _selective_merge_interval_lists(list_out, args)
        list_out.extend(_)
        for i in list_out:
            anno_counter += 1
            attr = f"ID={seqid}{strand_letter_func(strand)}_{anno_counter}" \
                   f";Name={seqid}_{strand_letter_func(strand)}_{new_type}_{anno_counter}" \
                   f";seq_len={i[1] - i[0] + 1}" \
                   f";condition={i[3]};best_mean_coverage={round(i[2], 2)}"
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


def _selective_merge_interval_lists(list_in, args):
    start = 0
    end = 1
    cmean = 2
    name = 3
    last = -1
    select_func = lambda x, y: 1 if x[cmean] > y[cmean] else 2 if y[cmean] > x[cmean] else 0
    list_in = sorted(list_in, key=lambda l: (l[0], l[1]))
    list_out = []
    overlap_indices = []
    for loc in list_in:
        if not list_out:
            # append first element
            list_out.append(loc)
            continue
        if loc[start] in range(list_out[last][start], list_out[last][end] + 1 + args.merge_range):
            if loc[end] - list_out[last][start] + 1 > args.max_len:
                # Don't merge overlap if the length will be exceeded
                selected = select_func(list_out[last], loc)
                if selected == 1:
                    pass
                elif selected == 2:
                    list_out[last] = loc
                else:
                    list_out.append(loc)
            else:
                # Safe to merge
                list_out[last][1] = max([loc[end], list_out[last][end]])
                list_out[last][2] = max([list_out[last][cmean], loc[cmean]])
                list_out[last][3] += f",{loc[name]}"
                overlap_indices.append(list_out.index(list_out[last]))
            continue
        # Check if the gap between non overlapping locations are too close
        if loc[start] - list_out[last][end] - 1 <= args.max_gap:
            checked = select_func(list_out[last], loc)
            if checked == 2:
                list_out[last] = loc
        else:
            list_out.append(loc)

    overlap_indices = list(set(overlap_indices))
    overlap_indices.sort()
    overlaps_list_out = [list_out[i] for i in overlap_indices]
    return overlaps_list_out, [i for i in list_out if i not in overlaps_list_out]


def convert_wiggle_obj_to_arr(wig_obj, args):
    step_size = int(args.min_len / 3)
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
    #merged_df = merged_df.round({merged_df.columns[-2]: 0, merged_df.columns[-1]: 0})
    #merged_df.loc[merged_df[merged_df.columns[-2]] <= args.ignore_coverage / 2, merged_df.columns[-2]] = 0.0
    #merged_df.loc[merged_df[merged_df.columns[-1]] <= args.ignore_coverage / 2, merged_df.columns[-1]] = 0.0
    for seqid in merged_df["variableStep_chrom"].unique():
        tmp_merged = merged_df[merged_df["variableStep_chrom"] == seqid].drop("variableStep_chrom", axis=1).copy()
        arr_dict[seqid] = np.absolute(tmp_merged.to_numpy(copy=True))
    return arr_dict, cond_name


def get_chrom_sizes(fasta_pathes):
    ret_list = []
    for fasta_path in fasta_pathes:
        print(f"==> Parsing reference sequence: {os.path.basename(fasta_path)}")
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            ret_list.append({"seqid": seq_record.id,
                             "size": len(seq_record.seq),
                             "fasta": os.path.basename(fasta_path)})
    return ret_list

def create_wiggle_obj(fpath, chrom_sizes):
    return Wiggle(fpath, chrom_sizes, is_len_extended=True)

main()