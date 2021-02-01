import argparse
import pandas as pd
from os import path
from Bio import SeqIO


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--refseq_in", required=True, help="", type=str)
parser.add_argument("--fasta_out", required=True, help="", type=str)
args = parser.parse_args()
fasta_parsed = SeqIO.parse(path.abspath(args.refseq_in), "fasta")
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
fasta_out_str = ""
f_seq = ""
r_seq = ""
for seq_record in fasta_parsed:
    f_seq = str(seq_record.seq)
    r_seq = str(seq_record.reverse_complement().seq)[::-1] # reverse of reverse complement
    gff_df = gff_df[(gff_df["type"] == "CDS") & (gff_df["seqid"] == seq_record.id)]
    for index, row in gff_df.iterrows():
        start = row['start'] - 1
        end = row['end'] - 1
        attr = parse_attributes(row['attributes'])
        if row['strand'] == "+":
            seq = f_seq[start:end]
            fasta_out_str += f">{attr['gene']};{attr['locus_tag']}\n{seq}\n"
        elif row['strand'] == "-":
            seq = r_seq[start:end]
            fasta_out_str += f">{attr['gene']};{attr['locus_tag']}\n{seq[::-1]}\n"
        else:
            print("Fatal error")
print("Writing fasta file...")
with open(path.abspath(args.fasta_out), "w") as f:
    f.write(f"{fasta_out_str}")
