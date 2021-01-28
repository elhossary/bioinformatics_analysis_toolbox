import pandas as pd
from os import path
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref_gff", required=True, help="", type=str)
    parser.add_argument("--replace_gff", required=True, help="", type=str)
    parser.add_argument("--locus_attr_name", default="locus_tag", help="", type=str)
    parser.add_argument("--gene_attr_name", default="name", help="", type=str)
    parser.add_argument("--out_gff", required=True, help="", type=str)
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    ref_gff_df = pd.read_csv(path.abspath(args.ref_gff), names=col_names, sep="\t", comment="#")
    ref_gff_df = ref_gff_df[ref_gff_df["type"] == "gene"]

    genes_dict = {}
    for indx in ref_gff_df.index:
        attr = parse_attributes(ref_gff_df.at[indx, 'attributes'])
        genes_dict[attr[args.locus_attr_name]] = attr[args.gene_attr_name]
    out_str = open(path.abspath(args.replace_gff), "r").read()

    for k, v in genes_dict.items():
        out_str = out_str.replace(k, v)

    with open(path.abspath(args.out_gff), "w") as f:
        f.write(out_str)



def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()