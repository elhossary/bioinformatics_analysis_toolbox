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
    parser.add_argument("--separated", action='store_true', help="", default=False)
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(args.gff_in, names=col_names, sep="\t", comment="#")
    gff_df.sort_values(["seqid", "start", "end"], inplace=True)
    fasta_parsed = SeqIO.parse(args.refseq_in, "fasta")
    str_dict = {}
    for record in fasta_parsed:
        for indx in gff_df.index:
            if gff_df.at[indx, "seqid"] == record.id:
                point = gff_df.at[indx, "start"] - 1
                if gff_df.at[indx, "strand"] == "+":
                    slice_start = point - args.range
                    slice_end = point
                    if slice_start < 0:
                        slice_start = 0
                elif gff_df.at[indx, "strand"] == "-":
                    slice_start = point
                    slice_end = point + args.range
                    if slice_end > len(record.seq):
                        slice_end = len(record.seq)
                else:
                    slice_start = None
                    slice_end = None
                    exit(1)
                seq_str = str(record.seq)
                if args.separated:
                    str_dict[f"{record.id}_{indx}"] = seq_str[slice_start: slice_end]
                else:
                    if record.id not in str_dict.keys():
                        str_dict[record.id] = ""
                    str_dict[record.id] += seq_str[slice_start: slice_end]
    out_str = ""
    for k in str_dict.keys():
        str_dict[k] = '\n'.join(textwrap.wrap(str_dict[k], 80))
    for k, v in str_dict.items():
        out_str += f">{k}\n{v}\n"
    with open(args.fasta_out, "w") as f:
        f.write(out_str)
main()