import argparse
import pandas as pd
import os
import glob
from wiggle_matrix import WiggleMatrix
from Bio import SeqIO
import numpy as np
from scipy.signal import find_peaks, peak_prominences
import multiprocessing as mp


def generate_step_factor_col(in_col, step_range, orientation):
    df = pd.DataFrame()
    df["scores"] = in_col
    df["mean_before"] = df["scores"]
    df["mean_after"] = df["scores"].shift(-(step_range + 1))
    df["mean_before"] = df["mean_before"].rolling(step_range).mean()
    df["mean_after"] = df["mean_after"].rolling(step_range).mean()
    if orientation == "f":
        df["step_factor"] = df["mean_before"] - df["mean_after"]
    elif orientation == "r":
        df["step_factor"] = df["mean_after"] - df["mean_before"]
    else:
        print("Error")
        exit(1)
    df[df["step_factor"] < 0] = 0.0
    df["step_factor"] = df["step_factor"].shift(1)
    df.fillna(0.0, inplace=True)
    df.rename(columns={"step_factor": in_col.name}, inplace=True)
    return df[in_col.name]

parser = argparse.ArgumentParser()
parser.add_argument("--refseq_in", required=True, help="", type=str)
parser.add_argument("--wiggle_files", required=True, help="", type=str, nargs="+")
parser.add_argument("--step_range", default=3, help="", type=int)
parser.add_argument("--gff_out", required=True, help="", type=str)
parser.add_argument("--processes", default="auto", help="", type=str)

args = parser.parse_args()
args.processes = str(6) #str(mp.cpu_count())

fasta_parsed = SeqIO.parse(os.path.abspath(args.refseq_in), "fasta")
seqid_list = [seq_record.id for seq_record in fasta_parsed]
wiggle_pathes = []
df_cols = ["seqid", "location"]
for item in args.wiggle_files:
    for sub_item in glob.glob(item):
        wiggle_pathes.append(os.path.abspath(sub_item))
wig_matrix = WiggleMatrix(args.refseq_in, wiggle_pathes, args.processes)
f_wiggles_matrix, r_wiggles_matrix = wig_matrix.get_matrix_by_orientation()
f_wiggles_matrix = f_wiggles_matrix[f_wiggles_matrix.seqid.isin(seqid_list)]
r_wiggles_matrix = r_wiggles_matrix[r_wiggles_matrix.seqid.isin(seqid_list)]


f_pool = mp.Pool(processes=8)
r_pool = mp.Pool(processes=8)
f_processes = []
r_processes = []

for seqid in seqid_list:
    for col in f_wiggles_matrix.columns:
        if col not in df_cols:
            f_processes.append(f_pool.apply_async(
                generate_step_factor_col,
                (f_wiggles_matrix[f_wiggles_matrix.seqid == seqid][col], args.step_range, "f", )))
    for col in r_wiggles_matrix.columns:
        if col not in df_cols:
            r_processes.append(r_pool.apply_async(
                generate_step_factor_col,
                (r_wiggles_matrix[r_wiggles_matrix.seqid == seqid][col].abs(), args.step_range, "r", )))
f_done_process = [p.get() for p in f_processes]
r_done_process = [p.get() for p in r_processes]
f_cols = []
r_cols = []
for col in f_done_process:
    f_wiggles_matrix[col.name] = col
    f_cols.append(col.name)
for col in r_done_process:
    r_wiggles_matrix[col.name] = col * -1
    r_cols.append(col.name)



#wig_matrix.write_matrix_to_wiggle_files(f_wiggles_matrix, "/home/muhoss/analysis_2020/Term-Seq_coverage/", "step_factor_coverage_")
#wig_matrix.write_matrix_to_wiggle_files(r_wiggles_matrix, "/home/muhoss/analysis_2020/Term-Seq_coverage/", "step_factor_coverage_")

f_wiggles_matrix["maxed_step_factor_forward"] = f_wiggles_matrix.loc[:, f_cols].max(axis=1)
r_wiggles_matrix["maxed_step_factor_reverse"] = r_wiggles_matrix.loc[:, r_cols].min(axis=1)
f_wiggles_matrix.drop(f_cols, axis=1, inplace=True)
r_wiggles_matrix.drop(r_cols, axis=1, inplace=True)
#wig_matrix.write_matrix_to_wiggle_files(f_wiggles_matrix, "/home/muhoss/analysis_2020/Term-Seq_coverage/", "")
#wig_matrix.write_matrix_to_wiggle_files(r_wiggles_matrix, "/home/muhoss/analysis_2020/Term-Seq_coverage/", "")
out_lst = []

for seqid in seqid_list:
    counter = 0
    f = f_wiggles_matrix[f_wiggles_matrix["seqid"] == seqid]["maxed_step_factor_forward"].tolist()
    f_peaks, _ = find_peaks(f, threshold=60)
    #f_peaks, _ = find_peaks(f, prominence=10)

    f_TTS = [i + 2 for i in f_peaks.tolist()]
    for i in f_TTS:
        counter += 1
        out_lst.append([seqid, "TSS_prediction", "TTS", i, i, ".", "+", ".",
                        f"Name=TTS_{seqid}_{counter}_F;ID=TTS_{seqid}_{counter}_F"])
    r = r_wiggles_matrix[r_wiggles_matrix["seqid"] == seqid]["maxed_step_factor_reverse"].abs().tolist()
    r_peaks, _ = find_peaks(r, threshold=60)
    #r_peaks, _ = find_peaks(r, prominence=10)
    r_TTS = [i - 1 for i in r_peaks.tolist()]
    for i in r_TTS:
        counter += 1
        out_lst.append([seqid, "TSS_prediction", "TTS", i, i, ".", "-", ".",
                        f"Name=TTS_{seqid}_{counter}_R;ID=TTS_{seqid}_{counter}_R"])

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
out_gff_df = pd.DataFrame(out_lst, columns=col_names)
out_gff_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)
print(out_gff_df.shape)


exit(0)
