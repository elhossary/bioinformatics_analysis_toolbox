import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--coord_tsv_in", required=True, help="", type=str)
parser.add_argument("--cds_specie", required=True, help="", type=str)
parser.add_argument("--gff_out", required=True, help="", type=str)
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
coord_cols = ["start", "end", "cds_start", "cds_end", "ref_len", "cds_len", "score", "seqid", "cds_name"]
coord_df = pd.read_csv(os.path.abspath(args.coord_tsv_in), sep="\t", comment="#", skiprows=4, names=coord_cols)
#coord_df.drop(["cds_start", "cds_end"], inplace=True)
coord_df["strand"] = ""
coord_df["source"] = "nucmer"
coord_df["type"] = "ortholog"
coord_df["phase"] = "."
coord_df.sort_values(["seqid", "start", "end", "score"], inplace=True)
x = 0
for indx in coord_df.index:
    x += 1
    coord_df.at[indx, "strand"] = "-" if coord_df.at[indx, "cds_start"] > coord_df.at[indx, "cds_end"] else "+"
    coord_df.at[indx, "attributes"] = \
        f"id={coord_df.at[indx, 'seqid']}_orth_{x}"\
        f";name={args.cds_specie}_ortholog_{coord_df.at[indx, 'cds_name']}" \
        f";cds_name={coord_df.at[indx, 'cds_name']}"
coord_df.drop(["cds_start", "cds_end", "ref_len", "cds_len", "cds_name"], inplace=True, axis=1)
coord_df = coord_df.reindex(columns=col_names)

coord_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)