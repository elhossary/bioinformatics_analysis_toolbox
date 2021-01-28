import pandas as pd
from os import path
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref_gff", required=True, help="", type=str)
    parser.add_argument("--replace_gff", required=True, help="", type=str)
    parser.add_argument("--ref_locus_attr_name", default="locus_tag", help="", type=str)
    parser.add_argument("--ref_gene_attr_name", default="name", help="", type=str)
    parser.add_argument("--replace_attr_name", required=True, help="", type=str)
    parser.add_argument("--out_gff", required=True, help="", type=str)
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    ref_gff_df = pd.read_csv(path.abspath(args.ref_gff), names=col_names, sep="\t", comment="#")
    ref_gff_df = ref_gff_df[ref_gff_df["type"] == "gene"]
    replace_gff_df = pd.read_csv(path.abspath(args.replace_gff), names=col_names, sep="\t", comment="#")

    genes_dict = {}
    for indx in ref_gff_df.index:
        attr = parse_attributes(ref_gff_df.at[indx, 'attributes'])
        genes_dict[attr[args.ref_locus_attr_name]] = attr[args.ref_gene_attr_name]

    for indx in replace_gff_df.index:
        attr = parse_attributes(replace_gff_df.at[indx, 'attributes'])
        add_str = ";associated_gene_name="
        for i in attr[args.replace_attr_name].split(","):
            if i in attr.keys():
                add_str += f"{genes_dict[i]},"

        replace_gff_df.at[indx, 'attributes'] += add_str[:-1]
    print("Writing GFF file...")
    replace_gff_df.to_csv(path.abspath(args.out_gff), sep="\t", header=False, index=False)



def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()