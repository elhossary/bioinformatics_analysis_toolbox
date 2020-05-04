import argparse
import pandas as pd
from os import path
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--refseq_in", required=True, help="", type=str)
parser.add_argument("--end_range", default=10, required=False, help="", type=int)
args = parser.parse_args()
fasta_parsed = SeqIO.parse(path.abspath(args.refseq_in), "fasta")
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
str_out = ""
print("Writing GFF file...")
f_seq = ""
r_seq = ""
t_count = 0
for index, row in gff_df.iterrows():
    for seq_record in fasta_parsed:
        if row['seqid'] == seq_record.id:
            f_seq = str(seq_record.seq)
            r_seq = str(seq_record.reverse_complement().seq)
            break
    if row['strand'] == "+":

        t_count = f_seq[int(row['end']) - args.end_range:int(row['end'])].count("T")
    elif row['strand'] == "-":
        t_count = r_seq[int(row['start']):int(row['start']) + args.end_range].count("T")
    else:
        print("Fatal error")
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
        f";seq_len={int(row['end']) - int(row['start']) + 1}" + \
        f";t_count={t_count}" + \
        "\n"
outfile = open(f"T_counted_{path.basename(args.gff_in)}", "w")
outfile.write(f"###gff-version 3\n{str_out}###")
outfile.close()