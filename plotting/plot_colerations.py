import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import math


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
    x_tick_size = math.ceil((max(df.cutoff.values.tolist()) - min(df.cutoff.values.tolist())) / 50)
    x_ticks = list(range(int(min(df.cutoff.values.tolist())), int(max(df.cutoff.values.tolist())) + x_tick_size, x_tick_size))
    y_ceil = calc_ceil(max(df.term_count.values.tolist()))
    y_tick_size = calc_ceil(int(y_ceil / 24))
    y_ticks = list(range(0, y_ceil + y_tick_size, y_tick_size))
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)
    plt.grid(True)
    plt.legend()
    fig.savefig(os.path.abspath(args.plot_path_out))


def calc_ceil(x):
    return int(math.ceil(x / 100.0)) * 100


main()
