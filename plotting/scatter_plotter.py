import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("--csv_in", required=True, help="", type=str)
#parser.add_argument("--data_columns", required=True, help="", type=str, nargs='+')
parser.add_argument("--x_column", required=True, help="", type=str)
parser.add_argument("--y_column", required=True, help="", type=str)
parser.add_argument("--classes_column", required=False, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
args = parser.parse_args()
csv_df = pd.read_csv(os.path.abspath(args.csv_in), sep="\t", comment="#")
all_sources = list(set(csv_df['source'].values.tolist()))
sns.set(rc={'figure.figsize':(16, 9)})
if args.classes_column is None:
    fig = sns.scatterplot(data=csv_df, x=args.x_column, y=args.y_column, hue='source',
                          palette='husl').set_title(f"{args.x_column} vs. {args.y_column} per class".replace("_", ""))
else:
    fig = sns.scatterplot(data=csv_df[csv_df['source'] == args.classes_column], x=args.x_column, y=args.y_column,
                          palette='husl').set_title(f"{args.x_column} vs. {args.y_column}".replace("_", ""))
fig.figure.savefig(f"{os.path.abspath(args.file_out)}")