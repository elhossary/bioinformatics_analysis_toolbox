import argparse
import pandas as pd
from os import path

def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
parser.add_argument("--type", required=True, help="", type=str, choices=["csv", "excel"])
args = parser.parse_args()

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
gff_df.reset_index(inplace=True)
for index, row in gff_df.iterrows():
    gff_df.at[index, 'attributes'] = parse_attributes(row['attributes'])
    gff_df.at[index, 'attributes']['df_index'] = index

attr_df = pd.DataFrame(gff_df.attributes.values.tolist())
attr_df.set_index("df_index", inplace=True)
gff_df = pd.merge(left=gff_df, right=attr_df, left_index=True, right_index=True).fillna("")
gff_df.drop(['index', 'attributes', 'score', 'phase'], axis=1, inplace=True)
if args.type == "csv":
    gff_df.to_csv(args.file_out, sep="\t", header=True, index=False)
