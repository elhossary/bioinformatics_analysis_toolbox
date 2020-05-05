import argparse
import pandas as pd
from os import path


parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--energy_values_in", required=True, help="", type=str)
parser.add_argument("--gff_out", required=True, help="", type=str)
args = parser.parse_args()

energy_values_file = open(path.abspath(args.energy_values_in), "r")
counter = 0
energy_values_list = []
tmp = ""
for line in energy_values_file.readlines():
    counter += 1
    tmp += line
    if counter == 3:
        energy_values_list.append(tmp.split("\n"))
        tmp = ""
        counter = 0
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
str_out = ""
for index, row in gff_df.iterrows():
    for item in energy_values_list:
        if row['seqid'] in item[0] and str(row['start']) in item[0] and \
                str(row['end']) in item[0] and row['strand'] in item[0]:
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
                f";energy_value={item[2].split(' ')[-1].replace('(', '').replace(')', '')}" + \
                "\n"
print("Writing GFF file...")
outfile = open(path.abspath(args.gff_out), "w")
outfile.write(f"{str_out}")
outfile.close()