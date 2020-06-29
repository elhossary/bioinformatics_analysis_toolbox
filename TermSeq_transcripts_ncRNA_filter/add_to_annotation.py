import pandas as pd
from os import path
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--offset", default=20, required=False, help="", type=int)
args = parser.parse_args()

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
str_out = ""
print("Writing GFF file...")
for index, row in gff_df.iterrows():
    if row['strand'] == "+":
        if row['start'] - args.offset > 0:
            row['start'] = row['start'] - args.offset
        else:
            row['start'] = 1
    elif row['strand'] == "-":
        row['end'] = row['end'] + args.offset
    else:
        print("fatal error")
        exit()
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
outfile = open(f"offset_added_{path.basename(args.gff_in)}", "w")
outfile.write(f"###gff-version 3\n{str_out}###")
outfile.close()