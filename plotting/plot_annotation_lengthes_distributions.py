import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--term_gff_in", required=True, help="", type=str)
parser.add_argument("--drna_gff_in", required=True, help="", type=str)
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

for path in [os.path.abspath(args.term_gff_in), os.path.abspath(args.drna_gff_in)]:
    gff_df = pd.read_csv(path, names=col_names, sep="\t", comment="#")
    gff_df["anno_len"] = gff_df["end"] - gff_df["start"] + 1
    bins = len(gff_df["anno_len"].unique().tolist())
    disc = gff_df["anno_len"].describe()
    values = gff_df["anno_len"].tolist()
    fig = plt.figure(figsize=(16, 9))
    plt.hist(values, bins=bins)
    plt.title(f"Annotations lengths distribution,\ncount: {disc['count']}, mean: {round(disc['mean'], 2)}, max: {disc['max']}, median: {disc['50%']}")
    plt.xlabel(f"Annotation length")
    plt.xticks(range(0, 301, 25))
    plt.ylabel("Frequency")
    plt.grid(True)
    fig.savefig(os.path.abspath(f"{path}.png"))

