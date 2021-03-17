import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--filter35", required=True, help="", type=str)
    parser.add_argument("--filter10", required=True, help="", type=str)
    parser.add_argument("--spacer", required=True, help="", type=int)
    parser.add_argument("--left_pad", required=True, help="", type=int)
    parser.add_argument("--right_pad", required=True, help="", type=int)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    filter_10_len = len(args.filter10)
    filter_35_len = len(args.filter35)
    window_size = filter_35_len + args.spacer + filter_10_len
    for indx in gff_df.index:
        gff_df.at[indx, "attributes"] = gff_df.at[indx, "attributes"].rstrip(";")
        attr = parse_attributes(gff_df.at[indx, "attributes"])
        seq = attr["sequence"]
        seq_Len = len(seq)
        window_start = 0
        add_str = ";motif_filter="
        res_lst = []
        loop_iter = seq_Len - window_size - args.right_pad

        for i in range(0, loop_iter):
            sub_seq_start = args.left_pad + i
            sub_seq_end = args.left_pad + i + window_size
            sub_seq = seq[sub_seq_start: sub_seq_end]
            box35_str = sub_seq[0: filter_35_len]
            box10_str = sub_seq[window_size - filter_10_len:]

            if box35_str == args.filter35 and box10_str == args.filter10:
                res_lst.append((True, True))
            elif box35_str != args.filter35 and box10_str == args.filter10:
                res_lst.append((False, True))
            elif box35_str == args.filter35 and box10_str != args.filter10:
                res_lst.append((True, False))
            else:
                res_lst.append((False, False))
            window_start += 1
        if (True, True) in res_lst:
            add_str += "both"
        else:
            if (False, True) in res_lst and (True, False) in res_lst:
                add_str += "box 10 OR box 35"
            else:
                if (False, True) in res_lst:
                    add_str += "box 10 only"
                elif (True, False) in res_lst:
                    print(76576)
                    add_str += "box 35 only"
                else:
                    add_str += "none"
        gff_df.at[indx, "attributes"] += add_str
    gff_df.to_csv(os.path.abspath(args.gff_in), sep="\t", header=False, index=False)


def parse_filter(str_in, spacer=0):
    dict_out = {}
    for i, c in enumerate(str_in):
        if c == "A" or c == "C" or c == "G" or c == "T":
            dict_out[i + spacer] = c

    return dict_out


def parse_attributes(attr_str):
    attr_str = attr_str.rstrip(";")
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


main()
