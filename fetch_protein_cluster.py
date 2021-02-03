import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv_in", required=True, help="", type=str)
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--attr_name", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    tsv_arr = np.genfromtxt(os.path.abspath(args.tsv_in), delimiter='\t', comments="#")
    print(tsv_arr)
main()