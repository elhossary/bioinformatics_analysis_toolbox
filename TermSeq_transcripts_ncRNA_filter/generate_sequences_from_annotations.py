import pandas as pd
from os import path
import argparse
from Bio import SeqIO
import textwrap


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


def wrap_text(str_in):
    str_out = ""
    for i in textwrap.wrap(str_in, 80):
        str_out += f"{i}\n"
    return str_out[:-1]

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--refseq_in", required=True, help="", type=str)
parser.add_argument("--dir_out", required=True, help="", type=str)
args = parser.parse_args()
fasta_parsed = SeqIO.parse(path.abspath(args.refseq_in), "fasta")
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")

for seq_record in fasta_parsed:
    f_seq = str(seq_record.seq)
    r_seq = str(seq_record.reverse_complement().seq)[::-1] # reverse of reverse complement
    str_out = ""
    for idx in gff_df[gff_df["seqid"] == seq_record.id].index:
        if gff_df.at[idx, "strand"] == "+":
            str_out += f">F_{gff_df.at[idx, 'seqid']}_{parse_attributes(gff_df.at[idx, 'attributes'])['name']}\n" \
                       f"{wrap_text(f_seq[gff_df.at[idx, 'start']: gff_df.at[idx, 'end']])}\n"
        elif gff_df.at[idx, "strand"] == "-":
            str_out += f">R_{gff_df.at[idx, 'seqid']}_{parse_attributes(gff_df.at[idx, 'attributes'])['name']}\n" \
                       f"{wrap_text(r_seq[gff_df.at[idx, 'start']: gff_df.at[idx, 'end']][::-1])}\n"
        else:
            print("fatal_error")
    outfile = open(f"{path.abspath(args.dir_out)}/{seq_record.id}.fasta", "w")
    outfile.write(f"{str_out}")
    outfile.close()

exit(0)
