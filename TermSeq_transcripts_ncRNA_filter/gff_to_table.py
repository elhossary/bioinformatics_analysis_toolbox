import argparse
import pandas as pd
from os import path
from sklearn import preprocessing

def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
parser.add_argument("--type", required=True, help="", type=str, choices=["csv", "excel"])
parser.add_argument("--scale_columns", required=False, help="", type=str, nargs='+')
parser.add_argument("--classes_column", required=False, help="", type=str)
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
scaled_columns = []
if args.scale_columns is not None:
    for column in args.scale_columns:
        scaled_columns.append(f"scaled_{column}")
    if args.classes_column is None:
        scaler = preprocessing.MinMaxScaler()
        scaled_df = pd.DataFrame(scaler.fit_transform(gff_df[args.scale_columns].astype(float).abs()),
                                 columns=scaled_columns)
        scaled_df["combined_scores"] = scaled_df.sum(axis=1)
        gff_df = pd.merge(left=gff_df, right=scaled_df, left_index=True, right_index=True).fillna("").round(2)
    else:
        scaler = preprocessing.MinMaxScaler()
        classes = list(set(gff_df['source'].values.tolist()))
        scaled_df = pd.DataFrame()
        for item in classes:
            tmp_df = pd.DataFrame(
                scaler.fit_transform(gff_df[gff_df['source'] == item][args.scale_columns].astype(float).abs()),
                columns=scaled_columns)

            scaled_df = scaled_df.append(tmp_df)
        scaled_df.reset_index(inplace=True)
        scaled_df["combined_scores"] = scaled_df.sum(axis=1)
        gff_df = pd.merge(left=gff_df, right=scaled_df, left_index=True, right_index=True).fillna("").round(2)

if args.type == "csv":
    gff_df.to_csv(args.file_out, sep="\t", header=True, index=False)
elif args.type == "excel":
    gff_df.to_excel(args.file_out, header=True, index=False)
