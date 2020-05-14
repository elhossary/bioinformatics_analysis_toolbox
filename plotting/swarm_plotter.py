import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("--csv_in", required=True, help="", type=str)
#parser.add_argument("--data_columns", required=True, help="", type=str, nargs='+')
parser.add_argument("--data_column", required=True, help="", type=str)
parser.add_argument("--classes_column", required=True, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
args = parser.parse_args()
csv_df = pd.read_csv(os.path.abspath(args.csv_in), sep="\t", comment="#")
sns.set(rc={'figure.figsize':(16, 9)})
fig = sns.swarmplot(data=csv_df, x=args.classes_column, y=args.data_column,
                    hue='source', palette='husl').set_title(f"{args.data_column} per class".replace("_", " "))
fig.figure.savefig(f"{os.path.abspath(args.file_out)}")