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
for item in args.wiggle_files:
    for sub_item in glob.glob(item):
        wiggle_pathes.append(os.path.abspath(sub_item))
f_wiggles_matrix, r_wiggles_matrix = WiggleMatrix(args.fasta_in, wiggle_pathes).get_matrix_by_orientation()
f_wiggles_cond = [col for col in f_wiggles_matrix.columns.tolist() if "seqid" != col != "location"]
r_wiggles_cond = [col for col in r_wiggles_matrix.columns.tolist() if "seqid" != col != "location"]
# minimize the pandas dataframe size for faster processing
f_slicing_lst = []
r_slicing_lst = []
args.step_range += 1
for i in gff_df[gff_df["strand"] == "+"]["start"].values.tolist():
    for step in range(1, args.step_range, 1):
        f_slicing_lst.append(i + step)
        f_slicing_lst.append(i - step)
    #f_slicing_lst.append(i)
f_slicing_lst.sort()
for i in gff_df[gff_df["strand"] == "-"]["end"].values.tolist():
    for step in range(1, args.step_range, 1):
        r_slicing_lst.append(i + step)
        r_slicing_lst.append(i - step)
    #r_slicing_lst.append(i)
r_slicing_lst.sort()
needed_seqid_list = gff_df["seqid"].unique()
# Slice unneeded seqids and locations
f_wiggles_matrix = f_wiggles_matrix[(f_wiggles_matrix["location"].isin(f_slicing_lst)) &
                                    (f_wiggles_matrix["seqid"].isin(needed_seqid_list))]
r_wiggles_matrix = r_wiggles_matrix[(r_wiggles_matrix["location"].isin(r_slicing_lst)) &
                                    (r_wiggles_matrix["seqid"].isin(needed_seqid_list))]
###########
# Generating conditions average column
f_wiggles_matrix["cond_mean"] = f_wiggles_matrix.loc[:, f_wiggles_cond].mean(axis=1)
r_wiggles_matrix["cond_mean"] = r_wiggles_matrix.loc[:, r_wiggles_cond].mean(axis=1)
# Converting the needed columns of matrix dataframe to numpy array for faster processing
f_wiggles_matrix = f_wiggles_matrix.loc[:, ["seqid", "location", "cond_mean"]].to_numpy()
r_wiggles_matrix = r_wiggles_matrix.loc[:, ["seqid", "location", "cond_mean"]].to_numpy()
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
        average_score_after = f_wiggles_matrix[(f_wiggles_matrix[:, 0] == tss_seqid) &
                                               (f_wiggles_matrix[:, 1] >= tss_loc + 1) &
                                               (f_wiggles_matrix[:, 1] <= tss_loc + args.step_range)][:, 2].mean()
        average_score_before = f_wiggles_matrix[(f_wiggles_matrix[:, 0] == tss_seqid) &
                                                (f_wiggles_matrix[:, 1] >= tss_loc - args.step_range) &
                                                (f_wiggles_matrix[:, 1] <= tss_loc - 1)][:, 2].mean()
        if average_score_before > 0 < average_score_after:
            step_factor = average_score_after / average_score_before
            step_height = average_score_after - average_score_before
        else:
            step_factor = 0
            step_height = 0
    elif tss_strand == "-":
        tss_loc = gff_df.at[i, "end"]
        average_score_before = abs(r_wiggles_matrix[(r_wiggles_matrix[:, 0] == tss_seqid) &
                                                    (r_wiggles_matrix[:, 1] >= tss_loc + 1) &
                                                    (r_wiggles_matrix[:, 1] <= tss_loc + args.step_range)][:, 2].mean())
        average_score_after = abs(r_wiggles_matrix[(r_wiggles_matrix[:, 0] == tss_seqid) &
                                                   (r_wiggles_matrix[:, 1] >= tss_loc - args.step_range) &
                                                   (r_wiggles_matrix[:, 1] <= tss_loc - 1)][:, 2].mean())
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
        f"{gff_df.at[i, 'attributes']};ave_step_height={step_height};ave_step_factor={step_factor}" + \
        "\n"
print("\nWriting GFF file...")
outfile = open(os.path.abspath(args.gff_out), "w")
outfile.write(f"{str_out}")
outfile.close()
exit(0)
