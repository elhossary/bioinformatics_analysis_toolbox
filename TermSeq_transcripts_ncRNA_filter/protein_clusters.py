import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys


def _expand_attributes_to_columns(expanded_df):

    for indx in expanded_df.index:
        attr_dict = {}
        comma_list = [item.split("=", maxsplit=1)
                      for item in expanded_df.at[indx, "attributes"].replace(";;", ";").split(";")]
        # TODO check for badly written attributes
        for item in comma_list:
            if len(item) == 2:
                attr_dict[item[0]] = item[1]
            elif len(item) == 1:
                print(f" Attribute without a value ignored")
            else:
                print(f" Some attributes was malformed and ignored near line {indx} - {item}")
                pass
        for k in attr_dict.keys():
            expanded_df.at[indx, k] = attr_dict[k]
    expanded_df.drop(["attributes"], axis=1, inplace=True)
    return expanded_df


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

parser = argparse.ArgumentParser()
parser.add_argument("--tsvs_in", required=True, help="", type=str, nargs="+")
parser.add_argument("--gff1_in", required=True, help="", type=str)
parser.add_argument("--gff2_in", required=True, help="", type=str)
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
tsvs_pathes = []
for item in args.tsvs_in:
    tsvs_pathes.extend(glob.glob(item))
gff1_df = pd.read_csv(os.path.abspath(args.gff1_in), names=col_names, sep="\t", comment="#")
gff2_df = pd.read_csv(os.path.abspath(args.gff2_in), names=col_names, sep="\t", comment="#")
gff1_df = gff1_df[gff1_df["attributes"].str.contains(";protein_id=")]
gff2_df = gff2_df[gff2_df["attributes"].str.contains(";protein_id=")]
gff1_df = _expand_attributes_to_columns(gff1_df)
gff2_df = _expand_attributes_to_columns(gff2_df)
product_df = pd.DataFrame(columns=["gff1_protein_id", "gff1_product", "gff2_protein_id", "gff2_product", "match"])

for indx in gff1_df.index:
    if gff1_df.at[indx, "product"] == np.nan:
        continue
    x_df = gff2_df[gff2_df["product"] == gff1_df.at[indx, "product"]]
    if x_df.empty or x_df.shape[0] > 1:
        if not x_df.empty:
            print(x_df.to_string())
        continue
    product_df = product_df.append({"gff1_protein_id": gff1_df.at[indx, "protein_id"],
                                    "gff1_product": gff1_df.at[indx, "product"],
                                    "gff2_protein_id": x_df.iat[0, x_df.columns.get_loc("protein_id")],
                                    "gff2_product": x_df.iat[0, x_df.columns.get_loc("product")],
                                    "match": "exact"}, ignore_index=True)
print(product_df.to_string())
exit()

tsv_col_names = ["cluster", "accession", "definition", "organism", "taxid", "length"]
tsvs_df = pd.DataFrame(columns=tsv_col_names)
for tsv_path in tsvs_pathes:
    x = pd.read_csv(os.path.abspath(tsv_path), names=tsv_col_names, sep="\t", comment="#")
    tsvs_df = tsvs_df.append(x)


tsvs_arr = tsvs_df.to_numpy()
clusters_df = pd.DataFrame(columns=tsv_col_names)
counter = 0
df_len = gff1_df.shape[0]
counter2 = 0
prot1_ids = [parse_attributes(gff1_df.at[indx, "attributes"])["protein_id"].split(".")[0] for indx in gff1_df.index]
prot2_ids = [parse_attributes(gff2_df.at[indx, "attributes"])["protein_id"].split(".")[0] for indx in gff2_df.index]
print(prot_ids)
x_arr = tsvs_arr[np.isin(tsvs_arr, prot_ids)]
print(x_arr)

