import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--plot_out", required=True, help="", type=str)
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")

gff_df["anno_len"] = gff_df["end"] - gff_df["start"] + 1
bins = len(gff_df["anno_len"].unique().tolist())
values = gff_df["anno_len"].tolist()
v_mean = gff_df["anno_len"].mean()
v_max = gff_df["anno_len"].max()
fig = plt.figure(figsize=(16, 9))
plt.hist(values, bins=bins)
plt.title(f"Annotations lengths distribution, average seq len: {round(v_mean, 2)}, max: {v_max}")
plt.xlabel(f"Annotation length")
plt.ylabel("Frequency")
plt.grid(True)
fig.savefig(os.path.abspath(args.plot_out))

