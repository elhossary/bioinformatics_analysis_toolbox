import argparse
import pandas as pd
from os import path
from Bio import SeqIO
from more_itertools import consecutive_groups


def get_longest_continuous_t(sequence_str):
    count = 0
    indices = [i for i, a in enumerate(sequence_str, 1) if a == "T"]
    indices.sort()
    signals = [list(group) for group in consecutive_groups(indices)]
    for i in signals:
        if len(i) > count:
            count = len(i)
    return count

parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--refseq_in", required=True, help="", type=str)
parser.add_argument("--end_range", default=10, required=False, help="", type=int)
parser.add_argument("--gff_out", required=True, help="", type=str)
args = parser.parse_args()
fasta_parsed = SeqIO.parse(path.abspath(args.refseq_in), "fasta")
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
str_out = ""

f_seq = ""
r_seq = ""
t_count = 0
longest_t_count = 0

for seq_record in fasta_parsed:
    for index, row in gff_df.iterrows():
        if row['seqid'] == seq_record.id:
            f_seq = str(seq_record.seq)
            r_seq = str(seq_record.reverse_complement().seq)
        else:
            continue
        if row['strand'] == "+":
            t_count = f_seq[int(row['end']) - args.end_range:int(row['end'])].count("T")
            longest_t_count = get_longest_continuous_t(f_seq[int(row['start']):int(row['end'])])
        elif row['strand'] == "-":
            t_count = r_seq[int(row['start']):int(row['start']) + args.end_range].count("T")
            longest_t_count = get_longest_continuous_t(r_seq[int(row['start']):int(row['end'])])
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
            f";last_{args.end_range}_bases_T_count={t_count}" + \
            f";longest_continuous_T={longest_t_count}" + \
            "\n"
print("Writing GFF file...")
outfile = open(path.abspath(args.gff_out), "w")
outfile.write(f"{str_out}")
outfile.close()