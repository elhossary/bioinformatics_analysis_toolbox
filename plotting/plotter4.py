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
csv_df = pd.read_csv(os.path.abspath(args.csv_in), sep="\t", comment="#")
all_sources = list(set(csv_df['source'].values.tolist()))
lists = {source:np.absolute(csv_df[csv_df['source'] == source][args.data_column].values) for source in all_sources}
x_axis = list(set(csv_df[args.data_column].values.tolist()))
bottom_lst = np.array(0 * len(x_axis))
fig = plt.figure(figsize=(16, 9))
for source in all_sources:
    heights = [list(lists[source]).count(i) for i in x_axis]
    plt.bar(x_axis, heights, label=source, bottom=bottom_lst)
    bottom_lst = np.array(heights)

plt.legend(prop={'size': 10})
plt.title(f"{args.data_column.replace('_', ' ')} distribution (stacked)")
plt.xlabel(f"Values")
plt.ylabel("Frequency")
plt.xticks(x_axis)
#plt.yticks(range(0, 7, 1))
plt.grid(True)
fig.savefig(f"{os.path.abspath(args.file_out)}")
