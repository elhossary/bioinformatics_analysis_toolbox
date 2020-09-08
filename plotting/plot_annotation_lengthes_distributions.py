import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--gff_files_in", required=True, help="", type=str)
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
pathes = []
for item in args.gff_files_in:
    for sub_item in glob.glob(item):
        pathes.append(os.path.abspath(sub_item))
for path in pathes:
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

