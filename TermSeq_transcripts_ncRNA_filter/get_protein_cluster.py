import pandas as pd
import argparse
import typing
import os
import glob
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table_file", required=True, type=str, help="GFF files (space separated)")
    parser.add_argument("--pc_file", required=True, type=str, help="")
    args = parser.parse_args()

    tbl_df = pd.read_csv(os.path.abspath(args.table_file), comment="#", sep="\t")
    tbl_df["accession"] = tbl_df["protein_id"]
    tbl_df["accession"].update(tbl_df["new_protein_id"])
    pc_cols = ["protein_cluster", "accession", "definition", "organism", "taxid", "length"]
    pc_df = pd.read_csv(os.path.abspath(args.pc_file), comment="#", sep="\t", names=pc_cols)
    ret_df = pd.merge(left=tbl_df, right=pc_df, on=["accession"], how='left')
    ret_df.fillna("", inplace=True)
    ret_df.sort_values(["protein_cluster"], inplace=True)
    export_file = os.path.dirname(args.table_file) + "/clustered_" + os.path.basename(args.table_file)
    ret_df.to_csv(export_file, sep='\t', index_label="index")

if __name__ == '__main__':
    main()