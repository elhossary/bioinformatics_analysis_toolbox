import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import glob


def main():
    # Params
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv_in", required=True, help="CSV file", type=str)
    parser.add_argument("--plot_path_out", required=True, help="Output file: Choose PNG or pdf", type=str)
    args = parser.parse_args()
    #
    df = pd.read_csv(os.path.abspath(args.csv_in), sep='\t').dropna()
    fig = plt.figure(figsize=(16, 9))
    plt.plot(df.cutoff.values.tolist(), df.term_count.values.tolist(), label='Terminators')
    plt.plot(df.cutoff.values.tolist(), df.srna_count.values.tolist(), label="sRNAs before merge")
    plt.plot(df.cutoff.values.tolist(), df.merged_srna_count.values.tolist(), label="sRNAs after merge")
    plt.title("Correlation between counts of terminators and sRNA at different cutoffs")
    plt.xlabel("Minimum Term-Seq coverage")
    plt.ylabel("Counts")
    plt.xticks(range(int(min(df.cutoff.values.tolist())) - 1, int(max(df.cutoff.values.tolist())) + 1, 1))
    plt.yticks(range(0, max(df.term_count.values.tolist()) + 200, 100))
    plt.grid(True)
    plt.legend()
    fig.savefig(os.path.abspath(args.plot_path_out))


main()
