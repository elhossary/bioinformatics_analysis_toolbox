#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
import os
import glob
from Bio import SeqIO
import re
import argparse

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
    fasta_seqids = {}
    for path in fasta_pathes:
        fasta_seqids[os.path.splitext(os.path.basename(path))[0]] = [rec.id for rec in SeqIO.parse(path, "fasta")]

    all_seqids = []
    for v in fasta_seqids.values():
        for i in v:
            if i not in all_seqids:
                all_seqids.append(i)
    for wig in wiggle_pathes:
        if args.output_dir is None:
            output_dir = os.path.dirname(wig)
        else:
            output_dir = args.output_dir
        print(f"Splitting file: {os.path.basename(wig)}")
        with open(wig, "r") as rf:
            header_text, content_dict = parse_wig_str(rf.read())
        for prefix, seqids in fasta_seqids.items():
            print(f"==> Writing {prefix}_{os.path.basename(wig)}")
            with open(f"{output_dir}/{prefix}_{os.path.basename(wig)}", "w") as wf:
                wf.write(f"{header_text}\n")
                for seqid in seqids:
                    for header in content_dict.keys():
                        if seqid in header:
                            wf.write(header)
                            wf.write(content_dict[header])
        if args.keep_missing:
            drop_headers = []
            for header in content_dict.keys():
                for seqid in all_seqids:
                    if seqid in header:
                        drop_headers.append(header)
            for item in drop_headers:
                del content_dict[item]
            if len(content_dict) > 0:
                print(f"==> Writing UNDEFINED_{os.path.basename(wig)}")
                with open(f"{output_dir}/UNDEFINED_{os.path.basename(wig)}", "w") as wf2:
                    for header in content_dict.keys():
                        wf2.write(header)
                        wf2.write(content_dict[header])


def parse_wig_str(in_str):
    ret_dict = {}
    header_text = in_str.split("\n", maxsplit=1)[0]
    in_str = in_str.replace(header_text + "\n", "")
    all_headers = re.findall(r'^.*chrom=.*$', in_str, flags=re.MULTILINE | re.IGNORECASE)
    splitters = ""
    for header in all_headers:
        splitters += header + "|"
    splitters = f"({splitters[:-1]})"
    split_str_list = re.split(rf"{splitters}", in_str, flags=re.MULTILINE | re.IGNORECASE)
    content_list = [i for i in split_str_list if i != '']
    for i in range(0, len(content_list), 2):
        ret_dict[content_list[i]] = content_list[i + 1]
    return header_text, ret_dict


main()
