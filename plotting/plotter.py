import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt

# Param
parser = argparse.ArgumentParser()
parser.add_argument("--csv_in", required=True, help="", type=str)
parser.add_argument("--data_column", required=False, help="", type=str)
parser.add_argument("--plot_type", required=False, help="", type=str)
parser.add_argument("--source", required=False, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
args = parser.parse_args()
csv_df = pd.read_csv(os.path.abspath(args.csv_in), sep="\t", comment="#")
fig = plt.figure(figsize=(16, 9))
if args.plot_type == "hist":
    if args.source != 'all':
        val_list = csv_df[csv_df['source'] == args.source][args.data_column].values.tolist()
    else:
        val_list = csv_df[args.data_column].values.tolist()
    bins = len(list(set(val_list)))
    plt.hist(val_list, bins=bins)
    plt.title(f"{args.data_column.replace('_', ' ')} distribution")
    plt.xlabel(f"Values")
    plt.ylabel("Frequency")
    plt.xticks(range(int(min(val_list)) - 1, int(max(val_list)) + 1, 1))
    #plt.yticks(range(0, 7, 1))
    plt.grid(True)
    fig.savefig(f"{os.path.abspath(args.file_out)}")
if args.plot_type == "bar":
    if args.source != 'all':
        x_axis = list(set(csv_df[csv_df['source'] == args.source][args.data_column].values.tolist()))
        heights = []
        for i in x_axis:
            heights.append(csv_df[csv_df['source'] == args.source][args.data_column].values.tolist().count(i))

    else:
        #all_sources = csv_df.source.unique().tolist()
        x_axis = list(set(csv_df[args.data_column].values.tolist()))

    heights = []
    for i in x_axis:
        heights.append(csv_df[args.data_column].values.tolist().count(i))

    plt.xticks(x_axis)
    step = 1
    if int(max(heights) / 20) < 1:
        step = 1
    else:
        step = int(max(heights) / 20)
    plt.yticks(range(0, max(heights), step))
    plt.bar(x_axis, heights)
    plt.title(f"{args.data_column.replace('_', ' ')} distribution")
    plt.xlabel(f"Values")
    plt.ylabel("Frequency")
    plt.grid(True)
    fig.savefig(f"{os.path.abspath(args.file_out)}")