from poly_t_stretch_finder import PolyTStretchFinder
import argparse
import os
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--fasta_in", required=True, help="", type=str)
parser.add_argument("--out_file", required=True, help="", type=str)
parser.add_argument("--min_len", required=True, help="", type=int)
parser.add_argument("--t_content", required=True, help="", type=float)
args = parser.parse_args()
obj = PolyTStretchFinder(args.fasta_in)
df = obj.find_stretches(args.min_len, args.t_content)
obj.write_to_gff(df, os.path.abspath(args.out_file))