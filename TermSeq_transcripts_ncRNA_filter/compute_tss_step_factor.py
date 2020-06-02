from statistics import mean
import argparse
import pandas as pd
import os
import glob
import sys
from wiggle_matrix import WiggleMatrix

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--fasta_in", required=True, help="", type=str)
parser.add_argument("--step_range", default=3, help="", type=int)
parser.add_argument("--gff_out", required=True, help="", type=str)
parser.add_argument("--wiggle_files", required=True, help="", type=str, nargs="+")
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
wiggle_pathes = []
str_out = ""
args.step_range += 1
for item in args.wiggle_files:
    for sub_item in glob.glob(item):
        wiggle_pathes.append(os.path.abspath(sub_item))
f_wiggles_matrix, r_wiggles_matrix = WiggleMatrix(args.fasta_in, wiggle_pathes).get_matrix_by_orientation()
f_wiggles_cond = [col for col in f_wiggles_matrix.columns.tolist() if "seqid" != col != "location"]
r_wiggles_cond = [col for col in r_wiggles_matrix.columns.tolist() if "seqid" != col != "location"]

f_slicing_lst = []
r_slicing_lst = []
for i in gff_df[gff_df["strand"] == "+"]["start"].values.tolist():
    for step in range(1, args.step_range, 1):
        f_slicing_lst.append(i + step)
        f_slicing_lst.append(i - step)
f_slicing_lst.sort()
for i in gff_df[gff_df["strand"] == "-"]["end"].values.tolist():
    for step in range(1, args.step_range, 1):
        r_slicing_lst.append(i + step)
        r_slicing_lst.append(i - step)
r_slicing_lst.sort()
f_wiggles_matrix = f_wiggles_matrix[f_wiggles_matrix["location"].isin(f_slicing_lst)]
r_wiggles_matrix = r_wiggles_matrix[r_wiggles_matrix["location"].isin(r_slicing_lst)]
gff_len = gff_df.shape[0]
for i in gff_df.index.tolist():
    sys.stdout.flush()
    sys.stdout.write("\r" + f"Computing TSS step height/factor progress: {round(i / gff_len * 100, 2)}%")
    tss_seqid = gff_df.at[i, "seqid"]
    tss_strand = gff_df.at[i, "strand"]
    tmp_lst_after = []
    tmp_lst_before = []
    step_factor = 0
    step_height = 0
    if tss_strand == "+":
        tss_loc = gff_df.at[i, "start"]
        for cond in f_wiggles_cond:
            tmp_lst_after.append(f_wiggles_matrix[(f_wiggles_matrix["seqid"] == tss_seqid) &
                                                  (f_wiggles_matrix["location"].between
                                                   (tss_loc + 1, tss_loc + args.step_range))][cond].mean())
            tmp_lst_before.append(f_wiggles_matrix[(f_wiggles_matrix["seqid"] == tss_seqid) &
                                                   (f_wiggles_matrix["location"].between
                                                    (tss_loc - args.step_range, tss_loc - 1))][cond].mean())
        average_score_after = mean(tmp_lst_after)
        average_score_before = mean(tmp_lst_before)
        if average_score_before > 0 < average_score_after:
            step_factor = average_score_after / average_score_before
            step_height = average_score_after - average_score_before
        else:
            step_factor = 0
            step_height = 0
    elif tss_strand == "-":
        tss_loc = gff_df.at[i, "end"]
        for cond in r_wiggles_cond:
            tmp_lst_before.append(r_wiggles_matrix[(r_wiggles_matrix["seqid"] == tss_seqid) &
                                                   (r_wiggles_matrix["location"].between
                                                   (tss_loc + 1, tss_loc + args.step_range))][cond].mean())
            tmp_lst_after.append(r_wiggles_matrix[(r_wiggles_matrix["seqid"] == tss_seqid) &
                                                  (r_wiggles_matrix["location"].between
                                                   (tss_loc - args.step_range, tss_loc))][cond].mean())
        average_score_after = mean([abs(v) for v in tmp_lst_after])
        average_score_before = mean([abs(v) for v in tmp_lst_before])
        if average_score_before > 0 < average_score_after:
            step_factor = average_score_after / average_score_before
            step_height = average_score_after - average_score_before
        else:
            step_factor = 0
            step_height = 0
    else:
        print("Fatal error: strand orientation problem")
        exit(1)
    str_out += \
        f"{gff_df.at[i, 'seqid']}\t" + \
        f"{gff_df.at[i, 'source']}\t" + \
        f"{gff_df.at[i, 'type']}\t" + \
        f"{int(gff_df.at[i, 'start'])}\t" + \
        f"{int(gff_df.at[i, 'end'])}\t" + \
        f"{gff_df.at[i, 'score']}\t" + \
        f"{gff_df.at[i, 'strand']}\t" + \
        f"{gff_df.at[i, 'phase']}\t" + \
        f"{gff_df.at[i, 'attributes']};ave_stepHeight={step_height};ave_stepFactor={step_factor}" + \
        "\n"
print("\nWriting GFF file...")
outfile = open(os.path.abspath(args.gff_out), "w")
outfile.write(f"{str_out}")
outfile.close()
exit(0)
