import argparse
import os
import pandas as pd


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


parser = argparse.ArgumentParser()
parser.add_argument("--gff1_in", required=True, help="", type=str)
parser.add_argument("--gff2_in", required=True, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff1_df = pd.read_csv(os.path.abspath(args.gff1_in), names=col_names, sep="\t", comment="#")
gff2_df = pd.read_csv(os.path.abspath(args.gff2_in), names=col_names, sep="\t", comment="#")

