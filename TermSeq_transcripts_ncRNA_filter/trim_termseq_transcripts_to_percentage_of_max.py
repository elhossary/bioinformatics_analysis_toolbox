import argparse
import pandas as pd
import os
from wiggletools.wiggle import Wiggle
from Bio import SeqIO
import sys
from more_itertools import consecutive_groups

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--f_wig_in", required=True, help="", type=str)
    parser.add_argument("--r_wig_in", required=True, help="", type=str)
    parser.add_argument("--nth", required=True, help="", type=int)
    parser.add_argument("--min_len", default=15, help="", type=int)
    parser.add_argument("--max_len", default=300, help="", type=int)
    parser.add_argument("--trim_ends_only", default=False, action="store_true")
    parser.add_argument("--keep_lower_slices", default=False, action="store_true")
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(os.path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    gff_df["start"] = gff_df["start"].astype(int)
    chrom_sizes = get_chrom_sizes([os.path.abspath(args.refseq_in)])
    f_wig_df = Wiggle(os.path.abspath(args.f_wig_in), chrom_sizes).get_wiggle()
    r_wig_df = Wiggle(os.path.abspath(args.r_wig_in), chrom_sizes).get_wiggle()
    remove_indecies = []
    gff_df_len = gff_df.shape[0]
    seqid_list = [x["seqid"] for x in chrom_sizes]
    invalid = 0
    slices_gff_df = pd.DataFrame(columns=col_names)
    all_seqid_slices_gff_df = pd.DataFrame(columns=col_names)
    for seqid in seqid_list:
        f_wig_df_slice = f_wig_df[f_wig_df["variableStep_chrom"] == seqid]
        r_wig_df_slice = r_wig_df[r_wig_df["variableStep_chrom"] == seqid]
        for idx in gff_df[gff_df["seqid"] == seqid].index:
            sys.stdout.flush()
            sys.stdout.write("\r" + f"Sequence ID {seqid} progress: {round(idx / gff_df_len * 100, 1)}%")
            tmp_df = None
            if gff_df.at[idx, "strand"] == "+":
                tmp_df = f_wig_df_slice[f_wig_df_slice["location"].between(
                    gff_df.at[idx, "start"], gff_df.at[idx, "end"])].copy()
            elif gff_df.at[idx, "strand"] == "-":
                tmp_df = r_wig_df_slice[r_wig_df_slice["location"].between(
                    gff_df.at[idx, "start"], gff_df.at[idx, "end"])].copy()
                tmp_df["score"] = tmp_df["score"].abs()
            else:
                print("Fatal error")
                exit(1)
            remove_limit = tmp_df["score"].max() * ((100 - args.nth) / 100)
            tmp_df = tmp_df[tmp_df["score"] > remove_limit]
            if tmp_df.empty:
                remove_indecies.append(idx)
                invalid += 1
            else:
                if args.trim_ends_only:
                    new_len = tmp_df["location"].max() - tmp_df["location"].min() + 1
                    if args.min_len <= new_len <= args.max_len:
                        gff_df.at[idx, "start"] = tmp_df["location"].min()
                        gff_df.at[idx, "end"] = tmp_df["location"].max()
                    else:
                        remove_indecies.append(idx)
                    continue
                    # ELSE trim_ends_only
                locs_list = tmp_df["location"].tolist()
                cons_locs_list = [list(group) for group in consecutive_groups(locs_list)]
                if len(cons_locs_list) == 1:
                    new_len = tmp_df["location"].max() - tmp_df["location"].min() + 1
                    if args.min_len <= new_len <= args.max_len:
                        gff_df.at[idx, "start"] = tmp_df["location"].min()
                        gff_df.at[idx, "end"] = tmp_df["location"].max()
                    else:
                        remove_indecies.append(idx)
                else:
                    cons_locs_list = [cons_loc for cons_loc in cons_locs_list
                                      if args.min_len <= max(cons_loc) - min(cons_loc) + 1 <= args.max_len]
                    if len(cons_locs_list) == 0:
                        remove_indecies.append(idx)
                    elif len(cons_locs_list) == 1:
                        gff_df.at[idx, "start"] = min(cons_locs_list[0])
                        gff_df.at[idx, "end"] = max(cons_locs_list[0])
                    else:
                        gff_df_row_copy = gff_df.iloc[idx].copy()
                        remove_indecies.append(idx)
                        attr = parse_attributes(gff_df_row_copy["attributes"])
                        anno_id = attr["id"]
                        anno_name = attr["name"]
                        counter = 0
                        slices_tmp_df = pd.DataFrame(columns=col_names)
                        for cons_loc in cons_locs_list:
                            counter += 1
                            gff_df_row_copy["start"] = min(cons_loc)
                            gff_df_row_copy["end"] = max(cons_loc)
                            cov_mean = tmp_df[tmp_df["location"].isin(cons_loc)]["score"].mean()
                            gff_df_row_copy["attributes"] =\
                                f"ID={anno_id}_slice_{counter};Name={anno_name}_slice_{counter}"
                            gff_df_row_copy["cov_mean"] = cov_mean
                            slices_tmp_df = slices_tmp_df.append(gff_df_row_copy, ignore_index=True)
                        if args.keep_lower_slices:
                            slices_gff_df = \
                                slices_gff_df.append(slices_tmp_df.drop("cov_mean", axis=1), ignore_index=True)
                        else:
                            slices_tmp_df.sort_values(by="cov_mean", ascending=False, inplace=True)
                            max_row = slices_tmp_df.iloc[0].drop("cov_mean")
                            slices_gff_df = slices_gff_df.append(max_row, ignore_index=True)
        all_seqid_slices_gff_df = all_seqid_slices_gff_df.append(slices_gff_df, ignore_index=True)
    print(f"\nTotal annotations removed: {len(remove_indecies)}\n\t"
          f"- Invalid: {invalid}\n\t- Too short/long: {len(remove_indecies) - invalid}")
    gff_df.drop(remove_indecies, axis=0, inplace=True)
    gff_df = gff_df.append(all_seqid_slices_gff_df, ignore_index=True)
    gff_df.to_csv(os.path.abspath(f"{args.gff_out}"), sep="\t", header=False, index=False)


def get_chrom_sizes(fasta_pathes):
    ret_list = []
    for fasta_path in fasta_pathes:
        print(f"==> Parsing reference sequence: {os.path.basename(fasta_path)}")
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            ret_list.append({"seqid": seq_record.id,
                             "size": len(seq_record.seq),
                             "fasta": os.path.basename(fasta_path)})
    return ret_list


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()
exit(0)
