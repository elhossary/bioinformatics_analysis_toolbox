import argparse
import pandas as pd
import numpy as np
import os
import glob
import sys
from Bio import SeqIO
from wiggletools.wiggle import Wiggle
from wiggletools.wiggle_matrix import WiggleMatrix



def get_chrom_sizes(fasta_pathes):
    ret_list = []
    for fasta_path in fasta_pathes:
        print(f"==> Parsing reference sequence: {os.path.basename(fasta_path)}")
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            ret_list.append({"seqid": seq_record.id,
                             "size": len(seq_record.seq),
                             "fasta": os.path.basename(fasta_path)})
    return ret_list

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--fasta_in", required=True, help="", type=str)
parser.add_argument("--step_range", default=3, help="", type=int)
parser.add_argument("--gff_out", required=True, help="", type=str)
parser.add_argument("--wiggle_files", required=True, help="", type=str, nargs="+")
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
f_gff_df = gff_df.loc[gff_df["strand"] == "+"].reset_index()
r_gff_df = gff_df.loc[gff_df["strand"] == "-"].reset_index()
wiggle_pathes = []

for item in args.wiggle_files:
    for sub_item in glob.glob(item):
        wiggle_pathes.append(os.path.abspath(sub_item))
chrom_sizes = get_chrom_sizes(os.path.abspath([args.fasta_in]))
parsed_wiggles = [Wiggle(wiggle_path, chrom_sizes).get_wiggle() for wiggle_path in wiggle_pathes]
f_wiggles_matrix, r_wiggles_matrix = WiggleMatrix(parsed_wiggles, chrom_sizes, processes=1).get_matrix_by_orientation()
f_wiggles_cond = [col for col in f_wiggles_matrix.columns.tolist() if "seqid" != col != "location"]
r_wiggles_cond = [col for col in r_wiggles_matrix.columns.tolist() if "seqid" != col != "location"]
# minimize the pandas dataframe size for faster processing

args.step_range += 1
needed_seqid_list = gff_df["seqid"].unique()
f_wiggles_matrix = f_wiggles_matrix.loc[f_wiggles_matrix.seqid.isin(needed_seqid_list)]
r_wiggles_matrix = r_wiggles_matrix.loc[r_wiggles_matrix.seqid.isin(needed_seqid_list)]
f_slicing_dict = {}
r_slicing_dict = {}
for seqid in needed_seqid_list:
    f_slicing_lst = []
    r_slicing_lst = []
    f_tss_locs = f_gff_df.loc[f_gff_df["seqid"] == seqid]["start"].values.tolist()
    r_tss_locs = r_gff_df.loc[r_gff_df["seqid"] == seqid]["end"].values.tolist()

    for i in f_tss_locs:
        for step in range(1, args.step_range, 1):
            f_slicing_lst.append(i + step)
            f_slicing_lst.append(i - step)
        #f_slicing_lst.append(i)
    f_slicing_lst.sort()
    f_slicing_dict[seqid] = f_slicing_lst
    for i in r_tss_locs:
        for step in range(1, args.step_range, 1):
            r_slicing_lst.append(i + step)
            r_slicing_lst.append(i - step)
        #r_slicing_lst.append(i)
    r_slicing_lst.sort()
    r_slicing_dict[seqid] = r_slicing_lst


# Slice unneeded seqids and locations
f_wiggles_matrix_sliced = pd.DataFrame()
r_wiggles_matrix_sliced = pd.DataFrame()
for k in f_slicing_dict.keys():
    tmp = f_wiggles_matrix.loc[(f_wiggles_matrix.location.isin(f_slicing_dict[k])) & (f_wiggles_matrix.seqid == k)]
    f_wiggles_matrix_sliced = f_wiggles_matrix_sliced.append(tmp, ignore_index=True)
    tmp = r_wiggles_matrix.loc[(r_wiggles_matrix.location.isin(r_slicing_dict[k])) & (r_wiggles_matrix.seqid == k)]
    r_wiggles_matrix_sliced = r_wiggles_matrix_sliced.append(tmp, ignore_index=True)
f_wiggles_matrix_sliced.reset_index(inplace=True)
r_wiggles_matrix_sliced.reset_index(inplace=True)
###########

f_wiggles_matrix_sliced["TSS_group"] = 0
r_wiggles_matrix_sliced["TSS_group"] = 0
f_wiggles_matrix_sliced["AB"] = ""
r_wiggles_matrix_sliced["AB"] = ""

f_gff_df_len = f_gff_df.shape[0]
print("Computing TSS step height/factor")
for i in f_gff_df.index:
    sys.stdout.flush()
    sys.stdout.write("\r" + f"==> for forward strand: {round(i / f_gff_df_len * 100)} %")
    tss_seqid = f_gff_df.at[i, "seqid"]
    tss_loc = f_gff_df.at[i, "start"]

    f_wiggles_matrix_sliced.loc[(f_wiggles_matrix_sliced["seqid"] == tss_seqid) &
                         (f_wiggles_matrix_sliced['location'].between(tss_loc + 1, tss_loc + args.step_range)),
                         ["TSS_group", "AB"]] = tss_loc, "A"
    f_wiggles_matrix_sliced.loc[(f_wiggles_matrix_sliced["seqid"] == tss_seqid) &
                         (f_wiggles_matrix_sliced['location'].between(tss_loc - args.step_range, tss_loc - 1)),
                         ["TSS_group", "AB"]] = tss_loc, "B"
print("\n")
r_gff_df_len = r_gff_df.shape[0]
for i in r_gff_df.index:
    sys.stdout.flush()
    sys.stdout.write("\r" + f"==> for reverse strand: {round(i / r_gff_df_len * 100)} %")
    tss_seqid = r_gff_df.at[i, "seqid"]
    tss_loc = r_gff_df.at[i, "end"]
    r_wiggles_matrix_sliced.loc[(r_wiggles_matrix_sliced["seqid"] == tss_seqid) &
                         (r_wiggles_matrix_sliced['location'].between(tss_loc + 1, tss_loc + args.step_range)),
                         ["TSS_group", "AB"]] = tss_loc, "B"
    r_wiggles_matrix_sliced.loc[(r_wiggles_matrix_sliced["seqid"] == tss_seqid) &
                         (r_wiggles_matrix_sliced['location'].between(tss_loc - args.step_range, tss_loc - 1)),
                         ["TSS_group", "AB"]] = tss_loc, "A"
print("\n")
f_wiggles_matrix_sliced.drop(["location"], axis=1, inplace=True)
r_wiggles_matrix_sliced.drop(["location"], axis=1, inplace=True)
f_wiggles_matrix_sliced = f_wiggles_matrix_sliced.groupby(["seqid", "TSS_group", "AB"], as_index=False).mean()
r_wiggles_matrix_sliced = r_wiggles_matrix_sliced.groupby(["seqid", "TSS_group", "AB"], as_index=False).mean()
f_wiggles_matrix_sliced["cond_mean"] = f_wiggles_matrix_sliced.loc[:, f_wiggles_cond].mean(axis=1)
r_wiggles_matrix_sliced["cond_mean"] = r_wiggles_matrix_sliced.loc[:, r_wiggles_cond].mean(axis=1)
f_wiggles_matrix_sliced.drop(f_wiggles_cond, axis=1, inplace=True)
r_wiggles_matrix_sliced.drop(r_wiggles_cond, axis=1, inplace=True)

f_wiggles_matrix_A = f_wiggles_matrix_sliced.loc[f_wiggles_matrix_sliced["AB"] == "A"]\
    .rename(columns={"TSS_group": "TSS_group_A", "cond_mean": "cond_mean_A"})
f_wiggles_matrix_A.drop(["AB"], axis=1, inplace=True)

f_wiggles_matrix_B = f_wiggles_matrix_sliced.loc[f_wiggles_matrix_sliced["AB"] == "B"]\
    .rename(columns={"TSS_group": "TSS_group_B", "cond_mean": "cond_mean_B"})
f_wiggles_matrix_B.drop(["AB"], axis=1, inplace=True)

r_wiggles_matrix_A = r_wiggles_matrix_sliced.loc[r_wiggles_matrix_sliced["AB"] == "A"]\
    .rename(columns={"TSS_group": "TSS_group_A", "cond_mean": "cond_mean_A"})
r_wiggles_matrix_A.drop(["AB"], axis=1, inplace=True)

r_wiggles_matrix_B = r_wiggles_matrix_sliced.loc[r_wiggles_matrix_sliced["AB"] == "B"]\
    .rename(columns={"TSS_group": "TSS_group_B", "cond_mean": "cond_mean_B"})
r_wiggles_matrix_B.drop(["AB"], axis=1, inplace=True)

f_wiggles_matrix_sliced = pd.merge(how='inner', left=f_wiggles_matrix_A, right=f_wiggles_matrix_B,
                                   left_on=["seqid", "TSS_group_A"], right_on=["seqid", "TSS_group_B"])
r_wiggles_matrix_sliced = pd.merge(how='inner', left=r_wiggles_matrix_A, right=r_wiggles_matrix_B,
                                   left_on=["seqid", "TSS_group_A"], right_on=["seqid", "TSS_group_B"])


f_wiggles_matrix_sliced.drop(["TSS_group_B"], axis=1, inplace=True)
r_wiggles_matrix_sliced.drop(["TSS_group_B"], axis=1, inplace=True)
f_wiggles_matrix_sliced.rename(columns={"TSS_group_A": "TSS"}, inplace=True)
r_wiggles_matrix_sliced.rename(columns={"TSS_group_A": "TSS"}, inplace=True)
f_wiggles_matrix_sliced.reset_index(inplace=True)
r_wiggles_matrix_sliced.reset_index(inplace=True)
# Calculating height / factor columns
f_wiggles_matrix_sliced["step_factor"] = \
    (f_wiggles_matrix_sliced["cond_mean_A"] / f_wiggles_matrix_sliced["cond_mean_B"])
f_wiggles_matrix_sliced["step_height"] = \
    (f_wiggles_matrix_sliced["cond_mean_A"] - f_wiggles_matrix_sliced["cond_mean_B"])
f_wiggles_matrix_sliced["step_factor"] = \
    f_wiggles_matrix_sliced["step_factor"].replace([np.inf, -np.inf, np.nan], 0).round(2).astype(str)
f_wiggles_matrix_sliced["step_height"] = \
    f_wiggles_matrix_sliced["step_height"].replace([np.inf, -np.inf, np.nan], 0).round(2).astype(str)

r_wiggles_matrix_sliced["step_factor"] = \
    (r_wiggles_matrix_sliced["cond_mean_A"].abs() / r_wiggles_matrix_sliced["cond_mean_B"].abs())
r_wiggles_matrix_sliced["step_height"] = \
    (r_wiggles_matrix_sliced["cond_mean_A"].abs() - r_wiggles_matrix_sliced["cond_mean_B"].abs())
r_wiggles_matrix_sliced["step_factor"] = \
    r_wiggles_matrix_sliced["step_factor"].replace([np.inf, -np.inf, np.nan], 0).round(2).astype(str)
r_wiggles_matrix_sliced["step_height"] = \
    r_wiggles_matrix_sliced["step_height"].replace([np.inf, -np.inf, np.nan], 0).round(2).astype(str)

f_gff_df = pd.merge(how='inner', left=f_gff_df, right=f_wiggles_matrix_sliced,
                    left_on=["seqid", "start"], right_on=["seqid", "TSS"])
r_gff_df = pd.merge(how='inner', left=r_gff_df, right=r_wiggles_matrix_sliced,
                    left_on=["seqid", "end"], right_on=["seqid", "TSS"])
out_df = f_gff_df.append(r_gff_df, ignore_index=True)
out_df.sort_values(["seqid", "start"], inplace=True)
out_df["attributes"] = out_df["attributes"] +\
                       ";ave_step_height=" + out_df["step_height"] + \
                       ";ave_step_factor=" + out_df["step_factor"]

out_df = out_df.loc[:, col_names]
out_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)
exit(0)