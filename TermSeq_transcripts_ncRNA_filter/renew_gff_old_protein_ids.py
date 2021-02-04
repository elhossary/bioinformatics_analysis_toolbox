import pandas as pd
from os import path
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--old_gff", required=True, help="", type=str)
    parser.add_argument("--new_gff", required=True, help="", type=str)
    parser.add_argument("--out_gff", required=True, help="", type=str)
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    old_gff_df = pd.read_csv(path.abspath(args.old_gff), names=col_names, sep="\t", comment="#")
    new_gff_df = pd.read_csv(path.abspath(args.new_gff), names=col_names, sep="\t", comment="#")
    new_gff_df = new_gff_df[(new_gff_df["attributes"].str.contains(";inference=")) & \
                            (new_gff_df["attributes"].str.contains(";protein_id="))].copy()
    ids_dict = {}
    names_dict = {}
    miss_counter = 0
    prod_dict = {}
    for indx in new_gff_df.index:
        attr = parse_attributes(new_gff_df.at[indx, "attributes"])
        new_size = new_gff_df.at[indx, "end"] - new_gff_df.at[indx, "start"]
        new_name = attr["gene"] if "gene" in attr.keys() else None
        ids_dict[get_old_prot_id(attr["inference"])] = [attr["protein_id"], new_size, new_name]
        if "gene" in attr.keys():
            names_dict[attr["gene"]] = [attr["protein_id"], new_size, attr["product"] if "product" in attr.keys() else None]
    for indx in old_gff_df.index:
        if ";protein_id=" not in old_gff_df.at[indx, "attributes"]:
            continue
        attr = parse_attributes(old_gff_df.at[indx, "attributes"])
        old_size = old_gff_df.at[indx, "end"] - old_gff_df.at[indx, "start"]
        old_name = attr["gene"]
        old_product = attr["product"]
        if attr["protein_id"] in ids_dict.keys():
            if ids_dict[attr["protein_id"]][1] == old_size and ids_dict[attr["protein_id"]][2] == old_name:
                old_gff_df.at[indx, "attributes"] += f";new_protein_id={ids_dict[attr['protein_id']][0]}"
        else:
            if "gene" in attr.keys():
                if attr["gene"] in names_dict.keys():
                    if names_dict[attr['gene']][1] == old_size or names_dict[attr['gene']][2] == old_product:
                        old_gff_df.at[indx, "attributes"] += f";new_protein_id={names_dict[attr['gene']][0]}"
                        continue
            miss_counter += 1
    print(f"{miss_counter} of {old_gff_df.shape[0]} protein IDs missed from translation")
    old_gff_df.to_csv(path.abspath(args.out_gff), sep="\t", header=False, index=False)


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


def get_old_prot_id(in_str):
    return in_str.split(":")[-1]

main()
