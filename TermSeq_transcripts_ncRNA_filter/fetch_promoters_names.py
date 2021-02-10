import pandas as pd
from os import path
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--attr_name", required=True, help="", type=str, nargs="+")
    parser.add_argument("--prom_tsv_in", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_in_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    prom_tsv_in_df =\
        pd.read_csv(path.abspath(args.prom_tsv_in), names=["Promoter", "Gene"], sep="\t", comment="#", skiprows=1)

    for x in gff_in_df.index:
        if any(word.lower() in gff_in_df.at[x, "attributes"].lower() for word in args.attr_name):
            attr = parse_attributes(gff_in_df.at[x, "attributes"])
            for y in prom_tsv_in_df.index:
                for i in args.attr_name:
                    attr_name = i.lower()
                    if attr_name in attr.keys() and attr[attr_name].lower() in prom_tsv_in_df.at[y, "Gene"].lower():
                        gff_in_df.at[x, "attributes"] += f";associated_{i}_promoter={prom_tsv_in_df.at[y, 'Promoter']}"
    gff_in_df.to_csv(path.abspath(args.gff_in), sep="\t", header=False, index=False)


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


main()
