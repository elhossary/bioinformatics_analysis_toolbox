import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--filter35", required=True, help="", type=str)
    parser.add_argument("--filter10", required=True, help="", type=str)
    parser.add_argument("--spacer", required=True, help="", type=int)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    filter_10 = parse_filter(args.filter10, args.spacer)
    filter_35 = parse_filter(args.filter35, 0)
    size_before = gff_df.shape[0]
    drop_lst = []
    for indx in gff_df.index:
        gff_df.at[indx, "attributes"] = gff_df.at[indx, "attributes"].rstrip(";")
        attr = parse_attributes(gff_df.at[indx, "attributes"])
        seq = attr["sequence"]
        flg10 = False
        flg35 = False
        for i, c in filter_10.items():
            if c == seq[i]:
                flg10 = True
            else:
                flg10 = False
                break
        for i, c in filter_35.items():
            if c == seq[i]:
                flg35 = True
            else:
                flg35 = False
                break
        add_str=";motif_="
        if flg10 and flg35:
            add_str += "both"
        elif not flg10 and flg35:
            add_str += "box 35 only"
        elif flg10 and not flg35:
            add_str += "box 10 only"
        else:
            add_str += "none"
        gff_df.at[indx, "attributes"] += add_str
    gff_df.to_csv(os.path.abspath(args.gff_in), sep="\t", header=False, index=False)


def parse_filter(str_in, spacer):
    dict_out = {}
    for i, c in enumerate(str_in):
        if c == "A" or c == "C" or c == "G" or c == "T":
            dict_out[i + spacer] = c

    return dict_out


def parse_attributes(attr_str):
    attr_str = attr_str.rstrip(";")
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


main()
