import argparse
import pandas as pd
import os
import glob
from wiggletools.wiggle import Wiggle
from wiggletools.wiggle_matrix import WiggleMatrix
from Bio import SeqIO
import numpy as np
from scipy.signal import find_peaks, peak_prominences
import multiprocessing as mp

parser = argparse.ArgumentParser()
parser.add_argument("--refseq_in", required=True, help="", type=str)
parser.add_argument("--wiggle_files", required=True, help="", type=str, nargs="+")
parser.add_argument("--gff_out", required=True, help="", type=str)
parser.add_argument("--processes", help="", type=int)
args = parser.parse_args()
fasta_parsed = SeqIO.parse(os.path.abspath(args.refseq_in), "fasta")
chrom_size = []
for seq_record in fasta_parsed:
    chrom_size.append({"seqid": seq_record.id, "size": len(seq_record.seq), "fasta": os.path.basename(args.refseq_in)})
wiggles_parsed = []
df_cols = ["seqid", "location"]
for item in args.wiggle_files:
    for sub_item in glob.glob(item):
        wiggles_parsed.append(Wiggle(os.path.abspath(sub_item), chrom_size).get_wiggle(is_full=False))
wig_matrix = WiggleMatrix(wiggles_parsed, chrom_size, processes=args.processes)
wig_matrix.get_matrix_by_orientation()
f_wiggles_matrix = wig_matrix.f_wiggle_matrix_df
r_wiggles_matrix = wig_matrix.r_wiggle_matrix_df

out_lst = []
seqid_list = wig_matrix.wiggle_matrix_df["seqid"].unique().tolist()
print("finding peaks")

for seqid in seqid_list:
    distance = 50
    height = 1
    prominence = 5
    counter = 0
    f = f_wiggles_matrix[f_wiggles_matrix["seqid"] == seqid]["agg_col_forward"].tolist()
    f_peaks, f_peaks_prop = find_peaks(f, distance=distance, prominence=prominence, height=height)
    f_TTS = [i + 1 for i in f_peaks.tolist()]
    for i in f_TTS:
        counter += 1
        out_lst.append([seqid, "TSS_prediction", "TTS", i, i, ".", "+", ".",
                        f"Name=TTS_{seqid}_{counter}_F;ID=TTS_{seqid}_{counter}_F"])
    r = r_wiggles_matrix[r_wiggles_matrix["seqid"] == seqid]["agg_col_reverse"].abs().tolist()
    r_peaks, r_peaks_prop = find_peaks(r, distance=distance, prominence=prominence, height=height)
    r_TTS = [i for i in r_peaks.tolist()]
    for i in r_TTS:
        counter += 1
        out_lst.append([seqid, "TSS_prediction", "TTS", i, i, ".", "-", ".",
                        f"Name=TTS_{seqid}_{counter}_R;ID=TTS_{seqid}_{counter}_R"])

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
out_gff_df = pd.DataFrame(out_lst, columns=col_names)
out_gff_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)
print(out_gff_df.shape)


exit(0)
