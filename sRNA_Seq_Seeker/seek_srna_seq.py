# Author: Muhammad Elhossary | elhossary@zbmed.de
import sys
from numpy import genfromtxt
import argparse
import glob
import os
import matplotlib.pyplot as plt
from gff_overlap_merger import GFF_Overlap_Merger as gff_mrg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tss_in", required=True, type=str, help="")
    parser.add_argument("--term_in", required=True, type=str, help="")
    parser.add_argument("--gff_out", required=True, type=str, help="")
    parser.add_argument("--max_len", required=True, type=int, help="")
    parser.add_argument("--min_len", required=True, type=int, help="")
    parser.add_argument("--merge_overlaps", action='store_true', help="")
    args = parser.parse_args()
    tss_arr = build_arr_form_gff(glob.glob(args.tss_in)[0])
    term_arr = build_arr_form_gff(glob.glob(args.term_in)[0])
    output_base_name = os.path.basename(args.gff_out)
    output_path = os.path.abspath(os.path.join(args.gff_out, os.pardir))
    print("\n\n--- sRNA Seq Seeker ---\n\n")
    print(f"Seeking for possible sRNA at sequences at length between {args.min_len} and {args.max_len}")
    srna_gff_str, term_matching_tss_counts, tss_matching_term_counts \
        = find_possible_sRNA(args.max_len, tss_arr, term_arr, args.min_len)
    plot_hist(term_matching_tss_counts, "How many TSSs are linked to each terminator",
              f"{output_path}/plot_TSS_to_Term_{output_base_name}.png")
    plot_hist(tss_matching_term_counts, "How many terminators are linked to each TSS",
              f"{output_path}/plot_Term_to_TSS_{output_base_name}.png")
    print("\nWriting output to file")
    outfile = open(args.gff_out, "w")
    outfile.write(f"###gff-version 3\n{srna_gff_str}###")
    outfile.close()
    if args.merge_overlaps:
        srna_gff_str, count_before, count_after = gff_mrg(srna_gff_str, 0, "sRNA").merge_overlaps()
        print("\nWriting merged output to file")
        print(f"Total annotations after merge: {count_after} of {count_before}")
        print(f"Merged ratio: {round((count_before - count_after) / count_before * 100, 2)}%")
        outfile = open(f"{output_path}/merged_{output_base_name}", "w")
        outfile.write(f"###gff-version 3\n{srna_gff_str}###")
        outfile.close())
    print("DONE")


def find_possible_sRNA(srna_max_length, tss_arr, term_arr, srna_min_length):
    # Number of Column names
    # 0 = accession
    # 1 = source
    # 2 = type
    # 3 = start
    # 4 = end
    # 5 = dot1
    # 6 = strand
    # 7 = dot2
    # 8 = attributes

    r_srna_gff_str = ""
    srna_count = 0
    tss_arr_len = len(tss_arr)
    tss_matching_term_counts = []
    term_matching_tss_counts = []
    for tss_index, tss_row in enumerate(tss_arr):
        tss_matching_term_counts.append(0)
        sys.stdout.flush()
        sys.stdout.write("\r" + f"Progress: {round(tss_index / tss_arr_len * 100, 2)}% | " +
                         f"{srna_count} possible sRNAs found      ...")
        for term_index, term_row in enumerate(term_arr):
            if tss_row[0] == term_row[0]:
                if tss_row[6] == term_row[6] == "+":

                    if tss_row[4] < term_row[3] and \
                            (term_row[4] - tss_row[3]) <= srna_max_length and \
                            srna_min_length <= (term_row[3] - tss_row[3]):
                        srna_count += 1
                        tss_matching_term_counts[-1] += 1
                        r_srna_gff_str += \
                            f"{tss_row[0]}\t" + \
                            f"sRNA_Seq_Seeker\t" + \
                            f"possible_sRNA_seq\t" + \
                            f"{tss_row[3]}\t" + \
                            f"{term_row[4]}\t" + \
                            f".\t" + \
                            f"{term_row[6]}\t" + \
                            f".\t" + \
                            f"id=possible_srna{srna_count};" + \
                            f"name=possible_srna{srna_count};" + \
                            f"seq_len={term_row[4] - tss_row[3]};" + \
                            f"matched_tss={parse_attributes(tss_row[8])['id']};" + \
                            f"matched_terminator={parse_attributes(term_row[8])['id']}\n"

                if tss_row[6] == term_row[6] == "-":

                    if term_row[4] < tss_row[3] and \
                            (tss_row[4] - term_row[3]) <= srna_max_length and \
                            srna_min_length <= (tss_row[3] - term_row[3]):
                        tss_matching_term_counts[-1] += 1
                        srna_count += 1
                        r_srna_gff_str += \
                            f"{tss_row[0]}\t" + \
                            f"sRNA_Seq_Seeker\t" + \
                            f"possible_sRNA_seq\t" + \
                            f"{term_row[3]}\t" + \
                            f"{tss_row[4]}\t" + \
                            f".\t" + \
                            f"{term_row[6]}\t" + \
                            f".\t" + \
                            f"id=possible_srna{srna_count};" + \
                            f"name=possible_srna{srna_count};" + \
                            f"seq_len={tss_row[4] - term_row[3]};" + \
                            f"matched_tss={parse_attributes(tss_row[8])['id']};" + \
                            f"matched_terminator={parse_attributes(term_row[8])['id']}\n"
    for term_index, term_row in enumerate(term_arr):
        term_matching_tss_counts.append(0)
        for tss_index, tss_row in enumerate(tss_arr):
            if tss_row[0] == term_row[0]:
                if tss_row[6] == term_row[6] == "+":
                    if tss_row[4] < term_row[3] and \
                            (term_row[4] - tss_row[3]) <= srna_max_length and \
                            srna_min_length <= (term_row[3] - tss_row[3]):
                        term_matching_tss_counts[-1] += 1
                if tss_row[6] == term_row[6] == "-":
                    if term_row[4] < tss_row[3] and \
                            (tss_row[4] - term_row[3]) <= srna_max_length and \
                            srna_min_length <= (tss_row[3] - term_row[3]):
                        term_matching_tss_counts[-1] += 1
    sys.stdout.write("\r" + f"Progress 100% with total {srna_count} possible sRNAs could be found")
    print("\n")
    return r_srna_gff_str, term_matching_tss_counts, tss_matching_term_counts


def plot_hist(list_in, title, output_file):
    distinct_list = []
    zero_counts = list_in.count(0)
    list_in = [i for i in list_in if i != 0]
    for i in list_in:
        if i not in distinct_list and i != 0:
            distinct_list.append(i)
    distinct_list.sort()
    bins = len(distinct_list)
    fig = plt.figure()
    plt.hist(list_in, bins=bins, rwidth=0.5)
    plt.title(f"{title}\nZero links count: {zero_counts}")
    plt.xlabel(f"Number of links")
    plt.ylabel("Frequency")
    plt.xticks(range(distinct_list[0], distinct_list[-1] + 1, 1))
    plt.grid(True)
    fig.savefig(output_file)


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


def build_arr_form_gff(path):
    data_arr = genfromtxt(path, delimiter="\t", comments="#", dtype=None, encoding=None)
    return data_arr

main()
