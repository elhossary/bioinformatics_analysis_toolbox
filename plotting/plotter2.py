import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--csv_in", required=True, help="", type=str)
parser.add_argument("--data_column", required=False, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
args = parser.parse_args()
csv_df = pd.read_csv(os.path.abspath(args.csv_in), sep="\t", comment="#").fillna(0)
all_sources = list(set(csv_df['source'].values.tolist()))
colors = ['red', 'green', 'blue']
lists = [np.absolute(csv_df[csv_df['source'] == source][args.data_column].values) for source in all_sources]
ticks = range(abs(int(csv_df[args.data_column].values.min() - 1)), abs(int(csv_df[args.data_column].values.max()) + 1), 1)
#bins = len(list(set(csv_df[args.data_column].values.tolist())))
bins = 50
fig = plt.figure(figsize=(16, 9))
plt.hist(lists, bins, stacked=True, density=False, label=all_sources)
plt.legend(prop={'size': 10})
plt.title(f"{args.data_column.replace('_', ' ')} distribution")
plt.xlabel(f"Values")
plt.ylabel("Frequency")
#plt.xticks(ticks)
#plt.yticks(range(0, 7, 1))
plt.grid(True)
fig.savefig(f"{os.path.abspath(args.file_out)}")
