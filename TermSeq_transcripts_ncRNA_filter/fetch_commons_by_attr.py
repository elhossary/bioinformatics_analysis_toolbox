import pandas as pd
from os import path
import argparse
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff1", required=True, help="", type=str)
    parser.add_argument("--gff2", required=True, help="", type=str)
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff1_df = pd.read_csv(path.abspath(args.gff1), names=col_names, sep="\t", comment="#")
    gff2_df = pd.read_csv(path.abspath(args.gff2), names=col_names, sep="\t", comment="#")

    gff2_df = gff2_df[gff2_df["attributes"].str.contains("ortholog")]
    gff1_df = gff1_df[gff1_df["attributes"].str.contains("ortholog")]
    names = []
    counter = 0
    for i in gff1_df.index:
        attr = parse_attributes(gff1_df.at[i, "attributes"])
        for k in attr.keys():
            if "downstream_cds_name" in k:
                names.append(attr[k].split(",")[0])

    names = set(names)
    orth_names = []

    for j in gff2_df.index:
        attr = parse_attributes(gff2_df.at[j, "attributes"])
        for k in attr.keys():
            if "_ortholog_name" in k:
                orth_names.append(attr[k].split(",")[0])
    orth_names = set(orth_names)
    print(len(names))
    print(len(orth_names))
    x = orth_names.intersection(names)
    print(len(x))

def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}
main()