#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
import os
import glob
from Bio import SeqIO
import re
import argparse
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq_files", required=True, help="", type=str, nargs="+")
    parser.add_argument("--wiggle_files", required=True, help="", type=str, nargs="+")
    parser.add_argument("--output_dir", default=None, help="", type=str)
    parser.add_argument("--keep_missing", action='store_true', help="")
    args = parser.parse_args()
    fasta_pathes = []
    for item in args.refseq_files:
        for sub_item in glob.glob(item):
            fasta_pathes.append(os.path.abspath(sub_item))
    wiggle_pathes = []
    for item in args.wiggle_files:
        for sub_item in glob.glob(item):
            wiggle_pathes.append(os.path.abspath(sub_item))
    fasta_seqids = [{os.path.splitext(os.path.basename(path))[0]:
                         list(set([rec.id for rec in SeqIO.parse(path, "fasta")]))} for path in fasta_pathes]

    unsplitted = ""
    for wig in wiggle_pathes:
        print(f"Splitting file: {os.path.basename(wig)}")
        with open(wig, "r") as rf:
            header_text, content_dict = parse_wig_str(rf.read())
            for fasta_seqids_list in fasta_seqids:
                for prefix, seqids in fasta_seqids_list.items():
                    if args.output_dir is None:
                        output_dir = os.path.dirname(wig)
                    else:
                        output_dir = args.output_dir
                    with open(f"{output_dir}/{prefix}_{os.path.basename(wig)}", "w") as wf:
                        wf.write(f"{header_text}\n")
                        for seqid in seqids:
                            for header in content_dict.keys():
                                if seqid in header:
                                    wf.write(header)
                                    wf.write(content_dict[header])
                                else:
                                    unsplitted += header
                                    unsplitted += content_dict[header]
                    if args.keep_missing:
                        with open(f"{output_dir}/UNDEFINED_{os.path.basename(wig)}", "w") as wf2:
                            wf2.write(f"{header_text}\n")
                            wf2.write(unsplitted)


def parse_wig_str(in_str):
    ret_dict = {}
    header_text = in_str.split("\n", maxsplit=1)[0]
    in_str = in_str.replace(header_text + "\n", "")
    all_headers = re.findall(r'^.*chrom=.*$', in_str, flags=re.MULTILINE | re.IGNORECASE)
    splitters = ""
    for header in all_headers:
        splitters += header + "|"
    splitters = f"{splitters[:-1]}"
    content_list = [i for i in re.split(rf"({splitters})", in_str, flags=re.MULTILINE | re.IGNORECASE) if i != '']
    for i in range(0, len(content_list), 2):
        ret_dict[content_list[i]] = content_list[i + 1]
    return header_text, ret_dict


main()
