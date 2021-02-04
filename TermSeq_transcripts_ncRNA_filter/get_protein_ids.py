import argparse
import os
import pandas as pd


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--file_out", required=True, help="", type=str)
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
gff_df = gff_df[(gff_df["attributes"].str.contains(";protein_id=")) | (gff_df["attributes"].str.contains(";new_protein_id="))]
prot_ids = []
for indx in gff_df.index:
    attr = parse_attributes(gff_df.at[indx, "attributes"])
    if "protein_id" in attr.keys():
        prot_ids.append(attr["protein_id"].split(".")[0])
    if "new_protein_id" in attr.keys():
        prot_ids.append(attr["new_protein_id"].split(".")[0])
prot_ids = list(set(prot_ids))
with open(os.path.abspath(args.file_out), "w") as f:
    f.write('\n'.join(prot_ids))
