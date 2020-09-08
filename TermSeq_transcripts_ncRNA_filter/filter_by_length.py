import pandas as pd
from os import path
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--min_len", default=25, required=False, help="", type=int)
parser.add_argument("--max_len", default=300, required=False, help="", type=int)
args = parser.parse_args()

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
abs_path = path.abspath(args.gff_in)
gff_df = pd.read_csv(abs_path, names=col_names, sep="\t", comment="#")
str_out = ""
print("Writing GFF file...")
for index, row in gff_df.iterrows():
    if args.min_len <= int(row['end']) - int(row['start']) + 1 <= args.max_len:
        str_out += \
            f"{row['seqid']}\t" + \
            f"{row['source']}\t" + \
            f"{row['type']}\t" + \
            f"{int(row['start'])}\t" + \
            f"{int(row['end'])}\t" + \
            f"{row['score']}\t" + \
            f"{row['strand']}\t" + \
            f"{row['phase']}\t" + \
            f"{row['attributes']}" + \
            "\n"
outfile = open(f"{path.dirname(abs_path)}/length_filtered_{path.basename(args.gff_in)}", "w")
outfile.write(f"###gff-version 3\n{str_out}###")
outfile.close()