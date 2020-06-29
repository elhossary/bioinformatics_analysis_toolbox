from Bio import SeqIO
from numpy import diff, where, split
import argparse
import pandas as pd
from wiggle_matrix import WiggleMatrix as WM
import glob

def main():
    # Params
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq_in", required=True, help="RefSeq fasta file", type=str)
    parser.add_argument("--tss_in", required=True, help="RefSeq fasta file", type=str)
    parser.add_argument("--wigs_in", required=True,
                        help="Term-Seq coverage file(s) (.wig), Must contain forward and reverse files", type=str)
    parser.add_argument("--gff_out", required=True, help="GFF output file name for terminators", type=str)
    parser.add_argument("--distance", required=True, help="Distance to look for terminator after a TSS", type=int)
    args = parser.parse_args()

    # ---------------------------
    print("Loading sequence file...")
    fasta_parsed = SeqIO.parse(glob.glob(args.refseq_in)[0], "fasta")
    wig_files = glob.glob(args.wigs_in)
    f_wigs_parsed, r_wigs_parsed = WM(wig_files, fasta_parsed).build_matrix()
    accession = ""

    # The following line is repeated due to the previous iterator exhaustion
    fasta_parsed = SeqIO.parse(glob.glob(args.refseq_in)[0], "fasta")
    for seq_record in fasta_parsed:
        f_seq_str = str(seq_record.seq)
        accession = seq_record.id
        print(f_wigs_parsed[accession].to_string())

def benchmark_terminators():
    pass
main()