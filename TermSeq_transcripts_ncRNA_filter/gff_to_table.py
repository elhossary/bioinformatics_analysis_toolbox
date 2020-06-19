import argparse
import pandas as pd
import numpy as np
from os import path
from sklearn import preprocessing


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


def add_excluded_values(df, exclude_list):
    exclude_list = ["scaled_" + item for item in exclude_list]
    exclude_list.append("combined_all_scores")
    cols = [col for col in df.columns.tolist() if col not in exclude_list]
    col_name = ""
    for c in cols:
        col_name += "_" + c.replace("scaled_", "")
    df[f"combined_score{col_name}"] = df[cols].sum(axis=1)
    return df

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
parser.add_argument("--type", required=True, help="", type=str, choices=["csv", "excel", "both"])
parser.add_argument("--scale_columns", required=False, help="", type=str, nargs='+')
parser.add_argument("--log_columns", required=False, help="", type=str, nargs='+')
parser.add_argument("--combine_exclude", required=False, help="", type=str, nargs='+')
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
if args.log_columns is not None:
    for log_col in args.log_columns:
        gff_df[f"{log_col}_no_log"] = gff_df[log_col]
        gff_df[log_col] = np.log10(gff_df[log_col].astype(float).replace([0, 0.0], np.nan))
scaled_columns = []
if args.scale_columns is not None:
    for column in args.scale_columns:
        scaled_columns.append(f"scaled_{column}")
    if args.classes_column is None:
        scaler = preprocessing.MinMaxScaler()
        scaled_df = pd.DataFrame(scaler.fit_transform(gff_df[args.scale_columns].to_numeric().astype(float).abs()),
                                 columns=scaled_columns)
        scaled_df["combined_all_scores"] = scaled_df.sum(axis=1)
        if args.combine_exclude is not None:
            scaled_df = add_excluded_values(scaled_df, args.combine_exclude)
        gff_df = pd.merge(left=gff_df, right=scaled_df, left_index=True, right_index=True).fillna("")
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
        scaled_df["combined_all_scores"] = scaled_df.sum(axis=1)
        if args.combine_exclude is not None:
            scaled_df = add_excluded_values(scaled_df, args.combine_exclude)
        gff_df = pd.merge(left=gff_df, right=scaled_df, left_index=True, right_index=True).fillna("")

gff_df = gff_df.round(2)

if args.type == "csv":
    gff_df.to_csv(path.abspath(f"{args.file_out}.csv"), sep="\t", header=True, index=False)
elif args.type == "excel":
    gff_df.to_excel(path.abspath(f"{args.file_out}.xlsx"), header=True, index=False)
elif args.type == "both":
    gff_df.to_csv(path.abspath(f"{args.file_out}.csv"), sep="\t", header=True, index=False)
    gff_df.to_excel(path.abspath(f"{args.file_out}.xlsx"), header=True, index=False)
