import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--target_gff_in", required=True, help="", type=str)
    parser.add_argument("--attr_id", required=True, help="", type=str)
    parser.add_argument("--target_attr_id", required=True, help="", type=str)
    parser.add_argument("--prefix", required=True, help="", type=str)
    parser.add_argument("--copy_target_attr", required=True, help="", type=str, nargs="+")
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), sep="\t", comment="#", names=col_names)
    target_gff_df = pd.read_csv(os.path.abspath(args.target_gff_in), sep="\t", comment="#", names=col_names)
    print(args.target_attr_id)
    print(args.attr_id)
    for indx in gff_df.index:
        if f";{args.attr_id}=" not in gff_df.at[indx, "attributes"]:
            continue
        row_attr = parse_attributes(gff_df.at[indx, "attributes"])
        find_word = f";{args.target_attr_id}={row_attr[args.attr_id]}"
        x_df = target_gff_df[target_gff_df["attributes"].str.contains(find_word)]
        if x_df.empty:
            continue
        x_dict = {}
        for x_indx in x_df.index:
            target_attr = parse_attributes(x_df.at[x_indx, "attributes"])
            for col in args.copy_target_attr:
                if col in target_attr.keys():
                    if f"{args.prefix}_{col}" not in x_dict.keys():
                        x_dict[f"{args.prefix}_{col}"] = []
                    x_dict[f"{args.prefix}_{col}"].append(target_attr[col])
        for k, v in x_dict.items():
            gff_df.at[indx, "attributes"] += f";{k}={','.join(v)}"
    gff_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)

def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()