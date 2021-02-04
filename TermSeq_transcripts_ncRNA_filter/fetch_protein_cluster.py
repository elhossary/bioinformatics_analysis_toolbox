import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv_in", required=True, help="", type=str)
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    tsv_col_names = ["cluster", "accession", "definition", "organism", "taxid", "length"]
    print("Loading TSV file...")
    tsv_arr = pd.read_csv(os.path.abspath(args.tsv_in), sep="\t", comment="#", names=tsv_col_names).to_numpy()
    print("TSV file Loaded")
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    gff_df_len = gff_df[gff_df["type"] == "CDS"].shape[0]
    counter = 0
    fetch_counter = 0
    for indx in gff_df.index:
        if gff_df.at[indx, "type"] != "CDS":
            continue
        counter += 1
        sys.stdout.flush()
        sys.stdout.write("\r" + f"Progress: {round(counter / gff_df_len * 100, 1)}% ")
        attr = parse_attributes(gff_df.at[indx, "attributes"])
        prot_id = attr["protein_id"] if "protein_id" in attr.keys() else None
        new_prot_id = attr["new_protein_id"] if "new_protein_id" in attr.keys() else None

        x_arr = tsv_arr[(tsv_arr[:, 1] == prot_id) | (tsv_arr[:, 1] == new_prot_id)]
        if x_arr.size > 0:
            gff_df.at[indx, "attributes"] += f";protein_cluster={x_arr[0, 0]}"
            fetch_counter += 1
    print(f"\n{fetch_counter} protein clusters could be fetched")
    gff_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)

def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()