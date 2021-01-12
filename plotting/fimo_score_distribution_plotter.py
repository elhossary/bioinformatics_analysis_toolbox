import pandas as pd
import argparse
import matplotlib.pyplot as plt
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fimo_tsv_file", required=True, help="", type=str)
    parser.add_argument("--out_file", required=True, help="", type=str)
    parser.add_argument("--title", help="", type=str)
    args = parser.parse_args()

    df = pd.read_csv(os.path.abspath(args.fimo_tsv_file), sep="\t")
    #df["score"] = df["score"].round(0)
    df["q-value"] = df["q-value"].astype(float)
    df["p-value"] = df["p-value"].astype(float)
    bins = len(df["q-value"].unique().tolist())
    disc = df["q-value"].describe()
    values = df["q-value"].tolist()
    fig = plt.figure(figsize=(16, 9))
    plt.hist(values, bins=bins)
    plt.title(f"{args.title}\n"
              f"count: {disc['count']}, mean: {round(disc['mean'], 2)}, max: {disc['max']}, median: {disc['50%']}")
    plt.xlabel(f"Q-value")
    plt.ylabel("Frequency")
    plt.grid(True)
    fig.savefig(os.path.abspath(f"{args.out_file}"))


main()
