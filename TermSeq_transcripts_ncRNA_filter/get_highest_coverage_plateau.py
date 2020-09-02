import argparse
import pandas as pd
import os
from wiggletools.wiggle import Wiggle
from Bio import SeqIO
import sys
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--f_wig_in", required=True, help="", type=str)
    parser.add_argument("--r_wig_in", required=True, help="", type=str)
    parser.add_argument("--plateau_size", required=True, help="", type=int)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    gff_df["start"] = gff_df["start"].astype(int)
    chrom_sizes = get_chrom_sizes([os.path.abspath(args.refseq_in)])
    f_wig_df = Wiggle(os.path.abspath(args.f_wig_in), chrom_sizes).get_wiggle()
    r_wig_df = Wiggle(os.path.abspath(args.r_wig_in), chrom_sizes).get_wiggle()
    gff_df_len = gff_df.shape[0]
    seqid_list = [x["seqid"] for x in chrom_sizes]
    plateaus_coverages = []
    remove_list = []
    for seqid in seqid_list:
        f_wig_df_slice = f_wig_df[f_wig_df["variableStep_chrom"] == seqid]
        r_wig_df_slice = r_wig_df[r_wig_df["variableStep_chrom"] == seqid]
        for idx in gff_df[gff_df["seqid"] == seqid].index:
            sys.stdout.flush()
            sys.stdout.write("\r" + f"Sequence ID {seqid} progress: {round(idx / gff_df_len * 100, 1)}%")
            anno_len = (gff_df.at[idx, "end"] - gff_df.at[idx, "start"] + 1)
            if anno_len < args.plateau_size:
                remove_list.append(idx)
                continue
            tmp_df = None
            if gff_df.at[idx, "strand"] == "+":
                tmp_df = f_wig_df_slice[f_wig_df_slice["location"].between(
                    gff_df.at[idx, "start"], gff_df.at[idx, "end"])].copy()
            elif gff_df.at[idx, "strand"] == "-":
                tmp_df = r_wig_df_slice[r_wig_df_slice["location"].between(
                    gff_df.at[idx, "start"], gff_df.at[idx, "end"])].copy()
                tmp_df["score"] = tmp_df["score"].abs()
            else:
                print("Fatal error")
                exit(1)
            max_score = tmp_df["score"].max()
            middle_max_loc = int(tmp_df[tmp_df['score'] == max_score]['location'].median())
            plateau_half_size = int(args.plateau_size / 2)
            mean_coverage = round(tmp_df[tmp_df["location"].between(
                middle_max_loc - plateau_half_size, middle_max_loc + plateau_half_size)]["score"].mean(), 2)
            gff_df.at[idx, "attributes"] += f";seq_len={anno_len};average_highest_plateau_coverage={mean_coverage}"\
                                            f";log10_AHPC={round(np.log10(mean_coverage), 2)}"
            plateaus_coverages.append(mean_coverage)
    gff_df.drop(remove_list, inplace=True)
    log_plateaus_coverages = [round(x, 2) for x in np.log10(plateaus_coverages)]
    unique_log_plateaus_coverages = []
    for i in log_plateaus_coverages:
        if i not in unique_log_plateaus_coverages:
            unique_log_plateaus_coverages.append(i)
    bins = len(unique_log_plateaus_coverages)
    fig = plt.figure(figsize=(16, 9))
    plt.hist(log_plateaus_coverages, bins=bins)
    plt.title(f"Distribution of highest plateau heights, plateau size: {args.plateau_size}")
    plt.xlabel(f"Average plateau coverage (rounded log10)")
    plt.xticks(range(0, int(max(unique_log_plateaus_coverages)) + 1, 1))
    plt.ylabel("Frequency")
    plt.grid(True)
    fig.savefig(f"{os.path.abspath(args.gff_out)}.png")
    gff_df.to_csv(os.path.abspath(f"{args.gff_out}"), sep="\t", header=False, index=False)

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
exit(0)