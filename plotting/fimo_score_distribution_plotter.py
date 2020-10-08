import pandas as pd
import argparse
import matplotlib.pyplot as plt
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fimo_tsv_file", required=True, help="", type=str)
    parser.add_argument("--out_file", required=True, help="", type=str)
    args = parser.parse_args()

    df = pd.read_csv(os.path.abspath(args.fimo_tsv_file), sep="\t")
    #df["score"] = df["score"].round(0)
    bins = len(df["score"].unique().tolist())
    disc = df["score"].describe()
    values = df["score"].tolist()
    fig = plt.figure(figsize=(16, 9))
    plt.hist(values, bins=bins)
    plt.title(f"FUR binding sites scores distribution\n"
              f"count: {disc['count']}, mean: {round(disc['mean'], 2)}, max: {disc['max']}, median: {disc['50%']}")
    plt.xlabel(f"FUR binding score")
    plt.ylabel("Frequency")
    plt.grid(True)
    fig.savefig(os.path.abspath(f"{args.out_file}"))


main()
