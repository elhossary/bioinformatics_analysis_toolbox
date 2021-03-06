import pandas as pd
from os import path
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    print("Appending old locus tags")
    for indx in gff_df.index:
        if "locus_tag" in gff_df.at[indx, "attributes"] and "old_locus_tag" not in gff_df.at[indx, "attributes"]:
            row_attr = parse_attributes(gff_df.at[indx, "attributes"])
            search_locus = row_attr["locus_tag"]
            x_df = gff_df[(gff_df["attributes"].str.contains(f"locus_tag={search_locus}")) &
                          (gff_df["attributes"].str.contains(f"old_locus_tag"))]
            if x_df.empty:
                continue
            for x_indx in x_df.index:
                old_locus = parse_attributes(x_df.at[x_indx, "attributes"])["old_locus_tag"]
                gff_df.at[indx, "attributes"] = gff_df.at[indx, "attributes"].strip(";")
                gff_df.at[indx, "attributes"] += f";old_locus_tag={old_locus}"
                break
    print("Writing GFF file")
    gff_df.to_csv(path.abspath(args.gff_in), sep="\t", header=False, index=False)


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()