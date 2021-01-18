import argparse
import pandas as pd
from Bio import SeqIO
import textwrap

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--range", default=50, help="", type=int)
    parser.add_argument("--fasta_out", required=True, help="", type=str)
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(args.gff_in, names=col_names, sep="\t", comment="#")
    gff_df.sort_values(["seqid", "start", "end"], inplace=True)
    fasta_parsed = SeqIO.parse(args.refseq_in, "fasta")
    str_dict = {}
    for record in fasta_parsed:
        str_dict[record.id] = ""

        for index in gff_df.index:
            if gff_df.at[index, "seqid"] == record.id:
                point = gff_df.at[index, "start"] - 1
                if gff_df.at[index, "strand"] == "+":
                    slice_start = point - args.range
                    slice_end = point
                    if slice_start < 0:
                        slice_start = 0
                elif gff_df.at[index, "strand"] == "-":
                    slice_start = point
                    slice_end = point + args.range
                    if slice_end > len(record.seq):
                        slice_end = len(record.seq)
                else:
                    slice_start = None
                    slice_end = None
                    exit(1)
                seq_str = str(record.seq)
                str_dict[record.id] += seq_str[slice_start: slice_end]
        str_dict[record.id] = '\n'.join(textwrap.wrap(str_dict[record.id], 80))
    out_str = ""
    print(str_dict)
    for k, v in str_dict.items():
        out_str += f">{k}\n{v}\n"
    with open(args.fasta_out, "w") as f:
        f.write(out_str)
main()