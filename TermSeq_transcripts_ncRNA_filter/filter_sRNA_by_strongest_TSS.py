import argparse
import pandas as pd
import os
from wiggletools.wiggle import Wiggle
from Bio import SeqIO
import numpy as np
import sys
from scipy import stats
import multiprocessing as mp
from more_itertools import consecutive_groups

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--f_tex_pos_wig_in", required=True, help="", type=str)
    parser.add_argument("--r_tex_pos_wig_in", required=True, help="", type=str)
    parser.add_argument("--f_tex_neg_wig_in", required=True, help="", type=str)
    parser.add_argument("--r_tex_neg_wig_in", required=True, help="", type=str)
    parser.add_argument("--nth", required=True, help="", type=int)
    parser.add_argument("--min_len", default=15, help="", type=int)
    parser.add_argument("--max_len", default=300, help="", type=int)
    parser.add_argument("--trim_ends_only", default=False, action="store_true")
    parser.add_argument("--keep_lower_slices", default=False, action="store_true")
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    gff_df["start"] = gff_df["start"].astype(int)
    gff_df["end"] = gff_df["end"].astype(int)
    chrom_sizes = get_chrom_sizes([os.path.abspath(args.refseq_in)])
    f_wig_df = Wiggle(os.path.abspath(args.f_wig_in), chrom_sizes).get_wiggle()
    r_wig_df = Wiggle(os.path.abspath(args.r_wig_in), chrom_sizes).get_wiggle()
    f_gff_df, r_gff_df = gff_df[gff_df["strand"] == "+"], gff_df[gff_df["strand"] == "-"]
    seqid_list = gff_df["seqid"].unique().tolist()
    grouped_gff_df = pd.DataFrame(columns=col_names + ["overlap_group"])
    for seqid in seqid_list:
        f_gff_df_seqid_slice = f_gff_df[f_gff_df["seqid"] == seqid].sort_values(by=["start", "end"])
        r_gff_df_seqid_slice = r_gff_df[r_gff_df["seqid"] == seqid].sort_values(by=["start", "end"])
        f_locs_list = [[f_gff_df_seqid_slice.at[x, "start"], f_gff_df_seqid_slice.at[x, "end"]]
                       for x in f_gff_df_seqid_slice.index]
        r_locs_list = [[r_gff_df_seqid_slice.at[x, "start"], r_gff_df_seqid_slice.at[x, "end"]]
                       for x in r_gff_df_seqid_slice.index]

        f_groups = merge_interval_lists(f_locs_list, 0, "all")
        r_groups = merge_interval_lists(r_locs_list, 0, "all")
        f_gff_df_seqid_slice["overlap_group"] = None
        r_gff_df_seqid_slice["overlap_group"] = None
        for x in f_gff_df_seqid_slice.index:
            for i in f_groups:
                if f_gff_df_seqid_slice.at[x, "start"] in range(i[0], i[1] + 1, 1):
                    f_gff_df_seqid_slice.at[x, "overlap_group"] = \
                        f"{r_gff_df_seqid_slice.at[x, 'seqid']}_{i[0]}..{i[1]}_+"
                    break
        for x in r_gff_df_seqid_slice.index:
            for i in r_groups:
                if r_gff_df_seqid_slice.at[x, "start"] in range(i[0], i[1] + 1, 1):
                    r_gff_df_seqid_slice.at[x, "overlap_group"] = \
                        f"{r_gff_df_seqid_slice.at[x, 'seqid']}_{i[0]}..{i[1]}_-"
                    break
        grouped_gff_df = grouped_gff_df.append(f_gff_df_seqid_slice)
        grouped_gff_df = grouped_gff_df.append(r_gff_df_seqid_slice)
        gff_df = filter_best_TSS(grouped_gff_df, f_wig_df, r_wig_df)
        gff_df.to_csv(os.path.abspath(f"{args.gff_out}"), sep="\t", header=False, index=False)


def filter_best_TSS(grouped_gff_df, f_wig_df, r_wig_df):
    groups_list = grouped_gff_df["overlap_group"].unique().tolist()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    out_df = pd.DataFrame(columns=col_names)
    for group in groups_list:
        tmp_df = grouped_gff_df[grouped_gff_df["overlap_group"] == group]
        tmp_df["TSS_height"] = None
        for x in tmp_df.index:
            if tmp_df.at[x, "strand"] == "+":

            elif tmp_df.at[x, "strand"] == "-":

            else:
                print("Fatal error in strand")
    return out_df

def get_chrom_sizes(fasta_pathes):
    ret_list = []
    for fasta_path in fasta_pathes:
        print(f"==> Parsing reference sequence: {os.path.basename(fasta_path)}")
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            ret_list.append({"seqid": seq_record.id,
                             "size": len(seq_record.seq),
                             "fasta": os.path.basename(fasta_path)})
    return ret_list

def merge_interval_lists(list_in, merge_range, annotate="all"):
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
    if annotate == 'all':
        return list_out
    overlap_indices = list(set(overlap_indices))
    overlap_indices.sort()
    if annotate == 'overlaps':
        list_out = [list_out[i] for i in overlap_indices]
    if annotate == 'no_overlaps':
        for i in reversed(overlap_indices):
            del list_out[i]
    return list_out


main()