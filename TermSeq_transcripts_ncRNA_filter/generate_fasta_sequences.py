import argparse
import pandas as pd
from os import path
from Bio import SeqIO


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--refseq_in", required=True, help="", type=str)
parser.add_argument("--end_range", default=40, required=False, help="", type=int)
parser.add_argument("--offset", default=5, required=False, help="", type=int)
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
    for index, row in gff_df.iterrows():
        if row['seqid'] == seq_record.id:
            if row['strand'] == "+":
                start = int(row['end']) - args.end_range - 1
                end = int(row['end']) + args.offset
                if start < int(row['start']):
                    start = int(row['start'])
                seq = f_seq[start:end].replace("T", "U")
                fasta_out_str += f">{row['seqid']}_{parse_attributes(row['attributes'])['name']}"\
                                 f":(+)_from_{start}_to_{end}\n{seq}\n"
            elif row['strand'] == "-":
                start = int(row['start']) - args.offset - 1
                end = int(row['start']) + args.end_range
                if int(row['end']) < end:
                    end = int(row['end'])
                seq = r_seq[start:end].replace("T", "U")
                fasta_out_str += f">{row['seqid']}_{parse_attributes(row['attributes'])['name']}"\
                                 f":(-)_from_{start}_to_{end}\n{seq[::-1]}\n"
            else:
                print("Fatal error")
print("Writing fasta file...")
outfile = open(path.abspath(args.fasta_out), "w")
outfile.write(f"{fasta_out_str}")
outfile.close()