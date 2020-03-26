import argparse
import glob
import os
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_in", required=True, type=str, help="")
    parser.add_argument("--f_wig_in", required=True, type=str, help="")
    #parser.add_argument("--r_wig_in", required=True, type=str, help="")
    args = parser.parse_args()
    fasta_parsed = SeqIO.parse(os.path.abspath(args.fasta_in), "fasta")
    seq_len = 0
    for seq_record in fasta_parsed:
        f_seq_str = str(seq_record.seq)
        accession = seq_record.id
        seq_len = len(seq_record.seq)


main()


