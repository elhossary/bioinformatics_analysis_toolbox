import pandas as pd
from os import path
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--fasta_in", required=True, help="", type=str)
parser.add_argument("--offset", default=20, required=False, help="", type=int)
parser.add_argument("--to", required=True, help="", type=str, choices=["start", "end", "both"])
args = parser.parse_args()

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
str_out = ""
f_gff_df = gff_df[gff_df["strand"] == "+"].copy()
r_gff_df = gff_df[gff_df["strand"] == "-"].copy()

if args.to == "start":
    f_gff_df.loc[:, ["start"]] = f_gff_df["start"] - args.offset
    r_gff_df.loc[:, ["end"]] = r_gff_df["end"] + args.offset
elif args.to == "end":
    r_gff_df.loc[:, ["start"]] = r_gff_df["start"] - args.offset
    f_gff_df.loc[:, ["end"]] = f_gff_df["end"] + args.offset
elif args.to == "both":
    f_gff_df.loc[:, ["start"]] = f_gff_df["start"] - args.offset
    f_gff_df.loc[:, ["end"]] = f_gff_df["end"] + args.offset
    r_gff_df.loc[:, ["end"]] = r_gff_df["end"] + args.offset
    r_gff_df.loc[:, ["start"]] = r_gff_df["start"] - args.offset
else:
    print("Fatal error")
    exit(0)
gff_df = f_gff_df.append(r_gff_df)
gff_df[gff_df["start"] < 0].loc[:, ["start"]] = 0
for seq_record in SeqIO.parse(path.abspath(args.fasta_in), "fasta"):
    genome_size = len(str(seq_record.seq))
    gff_df[(gff_df["seqid"] == seq_record.id) & (gff_df["end"] > genome_size)].loc[:, ["end"]] = genome_size
gff_df.sort_values(["seqid", "start", "end"], inplace=True)
print("Writing GFF file...")
out_file_name = f"{path.dirname(args.gff_in)}/offset_added_{path.basename(args.gff_in)}"
gff_df.to_csv(path.abspath(out_file_name), sep="\t", header=False, index=False)
