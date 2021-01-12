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
from statistics import mean
from more_itertools import consecutive_groups


def main():
    np.set_printoptions(suppress=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs="+")
    parser.add_argument("--merge_range", default=0, help="", type=int)
    parser.add_argument("--max_gap", default=20, help="", type=int)
    parser.add_argument("--peak_distance", default=40, help="", type=int)
    parser.add_argument("--min_len", default=50, help="", type=int)
    parser.add_argument("--max_len", default=350, help="", type=int)
    parser.add_argument("--ignore_coverage", default=10, help="", type=int)
    parser.add_argument("--sharpness", default=40, help="", type=int, choices=list(range(0, 101, 1)))
    parser.add_argument("--strict", default=False, help="", action='store_true')
    parser.add_argument("--step_size", default=5, help="", type=int)
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
    peaks_stats = peaks_gff_df["len"].describe(include='all')
    peaks_gff_df.drop("len", axis=1, inplace=True)
    print("Sorting annotated peaks")
    peaks_gff_df.sort_values(["seqid", "start", "end"], inplace=True)
    print(f"Total {peaks_stats['count']} peaks predicted, exporting to GFF")
          #f"median and average lengths: {peaks_stats['median']} and {peaks_stats['mean']}, exporting to GFF")
    peaks_gff_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)


def generate_locs(coverage_array, args, is_reversed, cond_name):
    location_col = 0
    raw_coverage_col = 1
    rising_col = 2
    falling_col = 3
    min_wid = (args.step_size * 2) - 1
    mean_func = lambda pos: mean(coverage_array[pos[0] - 1:pos[1] - 1, raw_coverage_col].tolist())
    rising_peaks, rising_peaks_props = find_peaks(coverage_array[:, rising_col], distance=args.min_len,
                                                  width=(min_wid, None), rel_height=1)
    falling_peaks, falling_peaks_props = find_peaks(coverage_array[:, falling_col], distance=args.min_len,
                                                    width=(min_wid, None), rel_height=1)
    rising_peaks_list = rising_peaks.tolist()
    falling_peaks_list = falling_peaks.tolist()
    falling_peaks_set = set(falling_peaks_list)
    rising_peaks_len = len(rising_peaks_list)
    counter = 0
    possible_locations = []
    for rp_id, rp in enumerate(rising_peaks_list):
        counter += 1
        sys.stdout.flush()
        sys.stdout.write("\r" + f"\t\tProgress: {round(counter / rising_peaks_len * 100, 1)}%")
        range_start = max(rp - args.max_len, 0) if is_reversed else rp + args.min_len
        range_end = rp - args.min_len if is_reversed else rp + args.max_len
        fp_range = set(range(range_start, range_end, 1))
        possible_fp = list(falling_peaks_set.intersection(fp_range))
        if not possible_fp:
            continue
        for fp in possible_fp:
            upper_loc = int(coverage_array[fp, location_col])
            lower_loc = int(coverage_array[rp, location_col])
            if is_reversed:
                upper_loc, lower_loc = lower_loc, upper_loc
            cmean = mean_func([lower_loc, upper_loc])
            if cmean == 0:
                continue
            lcratio = cmean / ((upper_loc - lower_loc) + 1)
            possible_locations.append([lower_loc, upper_loc, cmean, cond_name, lcratio])

    print("\n")
    # First round to filter the best annotation candidate for the same region (per condition)
    possible_locations = select_highest_mean_coverage_from_overlaps(possible_locations)
    # Second round to filter based on sharpness
    #possible_locations = filter_locations_based_on_sharpness(possible_locations, coverage_array, args, is_reversed)
    # Third round to connect annotations that are too close (gaps between fragments caused by read length bias)
    possible_locations = merge_locations_with_gaps(possible_locations, args, coverage_array)
    print(f"\t\tLocations found: {len(possible_locations)}")
    return possible_locations


def filter_locations_based_on_sharpness(list_in, coverage_array, args, is_reversed):
    print("\tFiltering based on required sharpness")
    if not list_in:
        return []
    list_out = []
    for loc in list_in:
        lower_loc = loc[0]
        upper_loc = loc[1]
        cmean = loc[2]
        lower_loc_id = lower_loc - 1
        upper_loc_id = upper_loc - 1
        shift = 1
        rp_bg_mean_cov = np.mean(coverage_array[lower_loc_id - (args.step_size + shift): lower_loc_id - shift, 1])
        fp_bg_mean_cov = np.mean(coverage_array[upper_loc_id + shift: upper_loc_id + shift + args.step_size, 1])
        if is_reversed:
            rp_bg_mean_cov, fp_bg_mean_cov = fp_bg_mean_cov, rp_bg_mean_cov
        rp_bg_cov_ratio = 100 - ((rp_bg_mean_cov / cmean) * 100)
        fp_bg_cov_ratio = 100 - ((fp_bg_mean_cov / cmean) * 100)
        if args.strict:
            if rp_bg_cov_ratio >= args.sharpness and fp_bg_cov_ratio >= args.sharpness:
                list_out.append(loc)
        else:
            if rp_bg_cov_ratio >= args.sharpness or fp_bg_cov_ratio >= args.sharpness:
                list_out.append(loc)
    return list_out


def select_highest_mean_coverage_from_overlaps(list_in):
    if not list_in:
        return []
    print("\tCleaning redundant locations")
    list_in = sorted(list_in, key=lambda l: (l[0], l[1]))
    start = 0
    end = 1
    cmean = 2
    lcratio = 4
    last = -1
    select_func = lambda x, y: 1 if x[cmean] > y[cmean] else 2 if y[cmean] > x[cmean] else 0

    list_out = [list_in[0]]
    for loc in list_in[1:]:
        if loc[start] == list_out[last][start] or loc[end] == list_out[last][end]:
            # If they share one of the ends, choose the the highest length coverage ratio
            if loc[lcratio] > list_out[last][lcratio]:
                list_out[last] = loc
            continue
        if loc[start] in range(list_out[last][start], list_out[last][end] + 2):
            # If they do not share any of the ends but overlapping
            selected = select_func(list_out[last], loc)
            if selected == 2 and loc[lcratio] > list_out[last][lcratio]:
                list_out[last] = loc
            else:
                continue
        else:
            list_out.append(loc)
    return list_out


def merge_locations_with_gaps(list_in, args, coverage_array):
    print("\tMerging locations with gaps")
    if not list_in:
        return []
    mean_func = lambda pos: mean(coverage_array[pos[0] - 1:pos[1] - 1, 1].tolist())
    max_mean_func = lambda x, y: 1 if x > y else y if x < y else 0
    gap_range = range(0, args.max_gap + 2, 1)
    length_range = range(args.min_len, args.max_len + 1, 1)
    list_in = sorted(list_in, key=lambda l: (l[0], l[1]))
    start = 0
    end = 1
    cmean = 2
    last = -1

    list_out = [list_in[0]]
    for loc in list_in[1:]:
        # check the gap between the new candidate and latest appended one
        if loc[start] - list_out[last][end] - 1 in gap_range:
            if loc[end] - list_out[last][start] + 1 in length_range:
                # If there is a gap and merged length is good
                if loc[cmean] / list_out[last][cmean] < 1 or list_out[last][cmean] / loc[cmean] < 1:
                    # This condition solves the issue of fragmented annotation because of the read biased coverage
                    # Checks if the average coverage between fragments is balanced
                    list_out[last][start] = min(list_out[last][start], loc[start])
                    list_out[last][end] = max(list_out[last][end], loc[end])
                    list_out[last][cmean] = mean_func([list_out[last][start], list_out[last][end]])
                else:
                    # if the average coverage between fragments is NOT balanced, choose one of them
                    max_mean = max_mean_func(list_out[last][cmean], loc[cmean])
                    if max_mean == 2:
                        list_out[last] = loc
            else:
                # if there is an allowed gap but the merged length is too long
                list_out.append(loc)
        else:
            # if no allowed gap
            list_out.append(loc)
    return list_out


def generate_annotations_from_positions(list_out, seqid, strand, new_type, args):
    strand_letter_func = lambda x: "F" if x == "+" else "R"
    anno_list = []
    if not list_out:
        return None
    list_out = merge_overlapping_locations(list_out, args)
    list_out = filter_locations_in_distance(list_out, args)
    anno_counter = 0
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


def merge_overlapping_locations(list_in, args):
    start = 0
    end = 1
    cmean = 2
    name = 3
    last = -1
    if not list_in:
        return []
    list_in = sorted(list_in, key=lambda l: (l[0], l[1]))
    list_out = [list_in[0]]
    for loc in list_in[1:]:
        if loc[start] in range(list_out[last][start], list_out[last][end] + 1) and \
                loc[end] - list_out[last][start] + 1 <= args.max_len:
            if list_out[last][end] - list_out[last][start] < loc[end] - loc[start]:
                list_out[last] = loc
            """
            list_out[last][end] = max([loc[end], list_out[last][end]])
            list_out[last][cmean] = max([list_out[last][cmean], loc[cmean]])
            list_out[last][name] += f",{loc[name]}"
            """
            continue
        list_out.append(loc)
    return list_out


def filter_locations_in_distance(list_in, args):
    print("\tFiltering unwanted locations")
    distance_range = range(0, args.peak_distance + 2, 1)
    list_in = sorted(list_in, key=lambda l: (l[0], l[1]))
    start = 0
    end = 1
    cmean = 2
    last = -1
    list_out = [list_in[0]]
    for loc in list_in[1:]:
        if loc[start] - list_out[last][end] - 1 in distance_range:
            if loc[cmean] > list_out[last][cmean]:
                list_out[last] = loc
                continue
        list_out.append(loc)
    return list_out


def convert_wiggle_obj_to_arr(wig_obj, args):
    arr_dict = {}
    wig_cols = ["variableStep_chrom", "location", "score"]
    wig_df = wig_obj.get_wiggle()
    cond_name = wig_df.iat[0, 1]
    wig_df = wig_df.loc[:, wig_cols]
    wig_df["score"] = wig_df["score"].abs()
    wig_df.loc[wig_df['score'] <= args.ignore_coverage, ["score"]] = 0.0
    merged_df = reduce(lambda x, y: pd.merge(x, y, on=["variableStep_chrom", "location"], how='left'),
                       [wig_df.loc[:, wig_cols],
                        wig_obj.to_step_height(args.step_size, "start_end").loc[:, wig_cols],
                        wig_obj.to_step_height(args.step_size, "end_start").loc[:, wig_cols]])
    merged_df["location"] = merged_df["location"].astype(int)

    for seqid in merged_df["variableStep_chrom"].unique():
        tmp_merged = merged_df[merged_df["variableStep_chrom"] == seqid].drop("variableStep_chrom", axis=1).copy()
        ret_arr = np.absolute(tmp_merged.to_numpy(copy=True))
        """
        for colid in [2, 3]:
            colname = tmp_merged.columns[colid]
            cgs = consecutive_groups(tmp_merged[tmp_merged[colname] > 0]["location"].tolist())
            for cg in cgs:
                cg_list = list(cg)
                start_id = min(cg_list) - 1
                end_id = max(cg_list) - 1
                sub_arr = ret_arr[start_id: end_id, colid]
                if sub_arr.size == 0:
                    continue
                target = np.where(ret_arr[start_id: end_id, colid] < np.max(sub_arr))
                ret_arr[target, colid] = 0.0
        """
        arr_dict[seqid] = ret_arr
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