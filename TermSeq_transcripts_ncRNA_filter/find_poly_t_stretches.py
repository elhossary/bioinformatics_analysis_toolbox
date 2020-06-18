from poly_t_stretch_finder import PolyTStretchFinder
import argparse
import os
import glob
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--fasta_in", required=True, help="", type=str)
parser.add_argument("--out_prefix", required=True, help="", type=str)
parser.add_argument("--min_len", required=True, help="", type=int)
parser.add_argument("--t_content", required=True, help="", type=int)
parser.add_argument()
args = parser.parse_args()
pathes = []

for item in args.fasta_in:
    for sub_item in glob.glob(item):
        pathes.append(os.path.abspath(sub_item))
obj = PolyTStretchFinder(pathes).find_stretches(args.min_len, args.t_content)
PolyTStretchFinder.write_to_gff(args.out_prefix, obj)
