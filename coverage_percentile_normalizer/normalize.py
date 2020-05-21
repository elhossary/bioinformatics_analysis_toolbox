from wiggle_parser import WiggleParser as wigprs
import argparse
import os
import numpy as np
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wigs_in", required=True, help="", type=str, nargs='+')
    parser.add_argument("--nth_percentile", required=True, help="", type=int)
    parser.add_argument("--dir_out", required=True, help="", type=str)
    args = parser.parse_args()
    wig_files_list = []
    for arg in [glob.glob(arg) for arg in args.wigs_in]:
        for item in arg:
            wig_files_list.append(os.path.abspath(item))

    for item in wig_files_list:
        wig_item_parsed = wigprs(item).parse()
        for key in wig_item_parsed:
            print(key)
            wig_percentile = get_percentile(wig_item_parsed[key], args.nth_percentile)
            print(wig_percentile)
            wig_item_parsed[key] = normalize(wig_item_parsed[key], wig_percentile)
        #wig_item_parsed(wig_item_parsed, args.dir_out)


def get_percentile(in_df, nth_percentile):
    return np.percentile(in_df[1].values.tolist(), nth_percentile)


def normalize(in_df, normalize_factor):
    in_df[1] = in_df[1] / normalize_factor
    return in_df


def write_file(in_dict, dir_path):
    for k in in_dict.keys():



main()