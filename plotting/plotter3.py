import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--csv_in", required=True, help="", type=str)
parser.add_argument("--data_column", required=False, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
parser.add_argument("--hist_type", required=True, help="", type=str, choices=["stacked", "overlaid"])
args = parser.parse_args()
csv_df = pd.read_csv(os.path.abspath(args.csv_in), sep="\t", comment="#")
all_sources = list(set(csv_df['source'].values.tolist()))
colors = ['black', 'orange', 'green']
lists = {source:np.absolute(csv_df[csv_df['source'] == source][args.data_column].values) for source in all_sources}
ticks = range(abs(int(csv_df[args.data_column].values.min() - 1)), abs(int(csv_df[args.data_column].values.max()) + 1), 1)
#bins = len(list(set(csv_df[args.data_column].values.tolist())))
#bins = 50
bins = np.linspace(0, 30, 50)
fig = plt.figure(figsize=(16, 9))

for index, k in enumerate(sorted(lists.keys())):
    if args.hist_type == "stacked":
        plt.hist(lists[k], bins, stacked=True, label=k, color=colors[index])
    elif args.hist_type == "overlaid":
        plt.hist(lists[k], bins, alpha=0.5, label=k, color=colors[index])
plt.legend(prop={'size': 10})
plt.title(f"{args.data_column.replace('_', ' ')} absolute distribution ({args.hist_type})")
plt.xlabel(f"Values")
plt.ylabel("Frequency")
#plt.xticks(range(0, 31, 1))
#plt.yticks(range(0, 7, 1))
plt.grid(True)
fig.savefig(f"{os.path.abspath(args.file_out)}")
