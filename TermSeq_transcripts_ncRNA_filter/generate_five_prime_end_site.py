import pandas as pd
from os import path
import argparse


def get_site(row):
    if row.strand == "+":
        row.end = row.start
    elif row.strand == "-":
        row.start = row.end
    else:
        print("FATAL ERROR")
        exit(1)
    return row


parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--gff_out", required=True, help="", type=str)
args = parser.parse_args()


col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
gff_df = gff_df.apply(lambda x: get_site(x), axis=1)

gff_df.to_csv(path.abspath(args.gff_out), sep="\t", index=False, header=False)
exit(0)
