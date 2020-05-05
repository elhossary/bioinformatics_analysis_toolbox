import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt

# Param
parser = argparse.ArgumentParser()
parser.add_argument("--csv_in", required=True, help="", type=str)
parser.add_argument("--data_column", required=False, help="", type=str)
args = parser.parse_args()

csv_df = pd.read_csv(os.path.abspath(args.csv_in), sep="\t", comment="#")
fig = plt.figure(figsize=(16, 9))
val_list = csv_df[args.data_column].values.tolist()
bins = len(list(set(val_list)))
plt.hist(val_list, bins=bins)
plt.title("Energy values distribution")
plt.xlabel(f"Values")
plt.ylabel("Frequency")
plt.xticks(range(int(min(val_list)) - 1, int(max(val_list)) + 1, 1))
plt.yticks(range(0, 25, 1))
plt.grid(True)
fig.savefig(f"{os.path.abspath(args.csv_in)}.png")