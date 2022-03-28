import pandas as pd
import argparse
from multiprocessing import Pool
import typing
import os
import glob
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gffs", required=True, type=str, nargs="+", help="GFF files (space separated)")
    args = parser.parse_args()
    parsed_paths = []
    for item in args.gffs:
        for sub_item in glob.glob(item):
            if not os.path.exists(os.path.abspath(sub_item)):
                print(f"Error: {sub_item} File does not exist!")
                sys.exit()
            parsed_paths.append(os.path.abspath(sub_item))
    for gff_path in parsed_paths:
        convert_to_csv(gff_path)
    sys.exit(0)


def convert_to_csv(in_path):
    column_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df = pd.read_csv(in_path, names=column_names, comment="#", sep="\t")
    if "attributes" not in df.columns:
        print("Warning: Attributes column not found!")
    else:
        df = explode_dict_yielding_func_into_columns(df, "attributes", parse_attributes)
        df.drop(["attributes"], inplace=True, axis=1)
    df.reset_index(inplace=True, drop=True)
    df.to_csv(f"{in_path}.tsv", sep='\t', index_label="index")


def explode_dict_yielding_func_into_columns\
                (df: pd.DataFrame, df_col: str, func: typing.Callable, column_prefix="") -> pd.DataFrame:
    seq_lst = df[df_col].tolist()
    df["TMP_COLUMN"] = pd.Series(list(map(func, seq_lst)))
    return explode_column_of_dicts(df, "TMP_COLUMN", column_prefix)


def explode_column_of_dicts(df: pd.DataFrame, df_col: str, column_prefix="") -> pd.DataFrame:
    df.reset_index(inplace=True, drop=True)
    tmp_df = df[df_col].apply(pd.Series)
    tmp_df.fillna("", inplace=True)
    if column_prefix != "":
        tmp_df = tmp_df.add_prefix(column_prefix)
    df = pd.merge(left=df, right=tmp_df, right_index=True, left_index=True, how='left')
    df.drop(columns=["TMP_COLUMN"], inplace=True)
    return df


def parse_attributes(attr_str: str):
    if attr_str == "":
        return {}
    attr_pairs = attr_str.split(";")
    attr_dict = {}
    for attr_pair in attr_pairs:
        attr_pair_lst = attr_pair.split("=")
        if len(attr_pair_lst) != 2:
            print(f"Warning: Skipping ambiguous key/value pair in GFF at: {attr_str}")
            continue
        k, v = attr_pair_lst[0], attr_pair_lst[1]
        if k.lower() in attr_dict.keys():
            if v in attr_dict[k.lower()]:
                continue
            attr_dict[k.lower()] += f"|{v}"
        else:
            attr_dict[k.lower()] = v
    return attr_dict


if __name__ == '__main__':
    main()
