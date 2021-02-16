import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--filter", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    filter_parsed = parse_filter(args.filter)
    drop_lst = []
    for indx in gff_df.index:
        try:
            attr = parse_attributes(gff_df.at[indx, "attributes"])
        except:
            continue
        seq = attr["sequence"]
        flg = False
        for i, c in filter_parsed.items():
            if c == seq[i]:
                flg = True
            else:
                flg = False
                break
        if not flg:
            drop_lst.append(indx)
    gff_df.drop(drop_lst, inplace=True, axis=0)
    gff_df.to_csv(os.path.abspath(args.gff_out), sep="\t", header=False, index=False)


def parse_filter(str_in):
    dict_out = {}
    for i, c in enumerate(str_in):
        if c == "A" or c == "C" or c == "G" or c == "T":
            dict_out[i] = c
    return dict_out


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


main()
