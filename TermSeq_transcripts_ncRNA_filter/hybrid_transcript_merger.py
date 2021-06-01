import io
import pandas as pd
import pybedtools as PyBed
import argparse
from itertools import product
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream_gff", required=True, type=str,
                        help="")
    parser.add_argument("--downstream_gff", required=True, type=str,
                        help="")
    parser.add_argument("--new_type_name", required=True, type=str,
                        help="")
    parser.add_argument("--window_size", default=400, type=int,
                        help="Maximum allowed annotation length")
    parser.add_argument("--prefix", default="", type=str,
                        help="")
    parser.add_argument("--gff_out", required=True, type=str,
                        help="Path to output GFF file")
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    upstream_df = pd.read_csv(os.path.abspath(args.upstream_gff), names=col_names, sep="\t", comment="#")
    downstream_df = pd.read_csv(os.path.abspath(args.downstream_gff), names=col_names, sep="\t", comment="#")
    seqid_list = [x for x in upstream_df["seqid"].unique() if x in downstream_df["seqid"].unique()]
    combinations = product(seqid_list, ["+", "-"])
    all_out_df = pd.DataFrame(columns=col_names)
    for comb in combinations:
        out_df = pd.DataFrame(columns=col_names)
        print(f"Processing {comb[0]}{comb[1]}")
        comb_upstream_df = \
            upstream_df[(upstream_df["seqid"] == comb[0]) & (upstream_df["strand"] == comb[1])].copy()
        comb_downstream_df = \
            downstream_df[(downstream_df["seqid"] == comb[0]) & (downstream_df["strand"] == comb[1])].copy()
        ret_dfs = merge_overlaps(comb_upstream_df, comb_downstream_df, args)
        out_df = out_df.append(ret_dfs[0], ignore_index=True)
        merged_imperfect, unmerged = merge_imperfect_overlaps(ret_dfs[1], comb[1], args)
        out_df = out_df.append(merged_imperfect, ignore_index=True)
        window_merged = merge_in_window(unmerged, args)
        out_df = subtract_overlaps(out_df, window_merged)
        out_df = merge_close_annotations(out_df, comb[1], args)
        all_out_df = all_out_df.append(out_df)
    all_out_df.sort_values(["seqid", "start", "end"], inplace=True)
    #out_df.drop_duplicates(subset=["seqid", "start", "end", "strand"], inplace=True)
    all_out_df.reset_index(inplace=True, drop=True)
    all_out_df = all_out_df.apply(lambda row: generate_attributes(row, args), axis=1)
    out_path = _inject_prefix(args.gff_out, args.prefix) if args.prefix != "" else os.path.abspath(args.gff_out)
    print(f"Exporting {all_out_df.shape[0]} annotations")
    all_out_df.to_csv(out_path, sep="\t", header=False, index=False)


def generate_attributes(row: pd.Series, args):
    strand_letter_func = lambda s: "F" if "+" in s else "R"
    row_counter = str(row.name)
    row["attributes"] = f"id={row['seqid']}_{strand_letter_func(row['strand'])}_{row_counter}" \
                        f";name={row['seqid']}_{strand_letter_func(row['strand'])}" \
                        f"_{args.new_type_name}_{row_counter}" \
                        f";seq_len={str(row['end'] - row['start'] + 1)}"
    return row


def merge_close_annotations(df, strand, args, dist=3):
    print("Merging close annotations")
    df_col = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    merge_col = ["seqid", "start", "end"]
    df = df[merge_col].copy()
    df_bed = PyBed.BedTool(df.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    merged_bed = df_bed.merge(d=dist)
    merged_df = pd.read_csv(io.StringIO(str(merged_bed)), names=merge_col, sep="\t", comment="#")
    merged_df["len"] = merged_df["end"] - merged_df["start"] + 1
    merged_df = merged_df[merged_df["len"] <= args.window_size]
    merged_df.drop(["len"], inplace=True, axis=1)
    merged_df["source"] = "ChimerAnno"
    merged_df["type"] = args.new_type_name
    merged_df["score"] = "."
    merged_df["strand"] = strand
    merged_df["phase"] = "."
    merged_df["attributes"] = ""
    return merged_df.reindex(columns=df_col)


def subtract_overlaps(df1, df2):
    df_col = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df1 = df1[df_col].copy()
    df2 = df2[df_col].copy()
    df1["attributes"] = "-"
    df2["attributes"] = "-"
    df1_bed = PyBed.BedTool(df1.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    df2_bed = PyBed.BedTool(df2.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    unique_bed = df2_bed.subtract(df1_bed, A=True)
    unique_df = pd.read_csv(io.StringIO(str(unique_bed)), names=df_col, sep="\t", comment="#")
    return df1.append(unique_df, ignore_index=True)


def merge_overlaps(up_df: pd.DataFrame, down_df: pd.DataFrame, args) -> tuple:
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "overlap_cluster"]
    merged_df = pd.DataFrame(columns=col_names)
    up_df["source"] = "upstream"
    down_df["source"] = "downstream"
    total_df = up_df.append(down_df, ignore_index=True)
    #up_df_bed = PyBed.BedTool(up_df.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    #down_df_bed = PyBed.BedTool(down_df.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    total_Bed = PyBed.BedTool(total_df.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    clustered_bed = total_Bed.cluster(d=-1)
    clustered_df = pd.read_csv(io.StringIO(str(clustered_bed)), names=col_names, sep="\t", comment="#")
    cluster_groups = clustered_df["overlap_cluster"].unique().tolist()
    print("Clustering overlaps")
    for cluster_group in cluster_groups:
        tmp_df = clustered_df[clustered_df["overlap_cluster"] == cluster_group].copy()
        clustered_df.drop(tmp_df.index, inplace=True)
        up_flg = True if "upstream" in tmp_df["source"].tolist() else False
        down_flg = True if "downstream" in tmp_df["source"].tolist() else False
        if up_flg and not down_flg:
            clustered_df = clustered_df.append(tmp_df, ignore_index=True)
            continue
        if not up_flg and down_flg:
            clustered_df = clustered_df.append(tmp_df, ignore_index=True)
            continue
        merge_return = _merge_perfect_overlaps_by_cluster(tmp_df, args)
        if type(merge_return) == pd.Series:
            merged_df = merged_df.append(merge_return, ignore_index=True)
        elif type(merge_return) == pd.DataFrame:
            clustered_df = clustered_df.append(merge_return, ignore_index=True)
        else:
            print("Internal error")
            exit(1)
    #unmerged_upstream_df.drop(["overlap_cluster"], inplace=True, axis=1)
    #unmerged_downstream_df.drop(["overlap_cluster"], inplace=True, axis=1)
    merged_df.drop(["overlap_cluster"], inplace=True, axis=1)
    return merged_df, clustered_df


def merge_in_window(df, args):
    print("Merging by window")
    merged_df = pd.DataFrame()
    df_cols = df.columns.tolist()
    #df_cols.remove("overlap_cluster")
    df = _extend_downstream_ends(df, args)
    gff_df = df[df_cols]
    gff_df_Bed = PyBed.BedTool(gff_df.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    print("Re-clustering in window")
    clustered_bed = gff_df_Bed.cluster(d=-1)
    clustered_df = pd.read_csv(io.StringIO(str(clustered_bed)), names=df_cols + ["window_cluster"], sep="\t", comment="#")
    clustered_df = clustered_df.merge(df, on=df_cols)
    up_clustered_df = clustered_df[clustered_df["source"] == "upstream"].copy()
    clustered_df.drop(up_clustered_df.index, inplace=True)
    up_clustered_df["start"] = up_clustered_df["original_start"]
    up_clustered_df["end"] = up_clustered_df["original_end"]
    clustered_df = clustered_df.append(up_clustered_df)
    clustered_df.drop(["original_start", "original_end"], inplace=True, axis=1)
    for cluster in clustered_df["window_cluster"].unique():
        tmp_df = clustered_df[clustered_df["window_cluster"] == cluster]
        clustered_df.drop(tmp_df.index, inplace=True)
        up_flg = True if "upstream" in tmp_df["source"].tolist() else False
        down_flg = True if "downstream" in tmp_df["source"].tolist() else False
        if (up_flg and not down_flg) or (not up_flg and down_flg):
            continue
        merged_df = merged_df.append(_merge_imperfect_window_overlaps(tmp_df, args), ignore_index=True)
    return merged_df


def merge_imperfect_overlaps(df: pd.DataFrame, strand, args):
    df_col = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    merge_col = ["seqid", "start", "end"]
    df = df[df_col]
    up_df = df[df["source"] == "upstream"]
    down_df = df[df["source"] == "downstream"]
    # BEDs
    up_down_bed = PyBed.BedTool(df.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    up_df_bed = PyBed.BedTool(up_df.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    down_df_bed = PyBed.BedTool(down_df.to_csv(sep="\t", header=False, index=False), from_string=True).sort()
    merged_up_down_bed = up_down_bed.merge(d=-1)
    fixed_merged_up_down_bed = merged_up_down_bed.intersect(up_df_bed)  # Trim merges based on upstream
    up_wo_down_bed = up_df_bed.subtract(down_df_bed, A=True)
    down_wo_up_bed = down_df_bed.subtract(up_df_bed, A=True)
    fixed_merged_up_down_bed = fixed_merged_up_down_bed.subtract(up_wo_down_bed, A=True)
    fixed_merged_up_down_bed = fixed_merged_up_down_bed.subtract(down_wo_up_bed, A=True)
    unmerged_up_bed = up_df_bed.intersect(fixed_merged_up_down_bed, v=True)
    unmerged_down_bed = down_df_bed.intersect(fixed_merged_up_down_bed, v=True)
    # DFs
    fixed_merged_up_down_df = \
        pd.read_csv(io.StringIO(str(fixed_merged_up_down_bed)), names=merge_col, sep="\t", comment="#")
    fixed_merged_up_down_df["source"] = "ChimerAnno"
    fixed_merged_up_down_df["type"] = args.new_type_name
    fixed_merged_up_down_df["score"] = "."
    fixed_merged_up_down_df["strand"] = strand
    fixed_merged_up_down_df["phase"] = "."
    fixed_merged_up_down_df["attributes"] = ""
    fixed_merged_up_down_df = fixed_merged_up_down_df.reindex(columns=df_col)
    unmerged_up_df = pd.read_csv(io.StringIO(str(unmerged_up_bed)), names=df_col, sep="\t", comment="#")
    unmerged_down_df = pd.read_csv(io.StringIO(str(unmerged_down_bed)), names=df_col, sep="\t", comment="#")
    unmerged_df = unmerged_up_df.append(unmerged_down_df, ignore_index=True)
    return fixed_merged_up_down_df, unmerged_df


def _merge_imperfect_window_overlaps(df: pd.DataFrame, args) -> pd.DataFrame:
    df = df.copy()
    is_reversed = True if "-" in df["strand"].tolist() else False
    ret_df = pd.DataFrame(columns=df.columns)
    ret_merge = df.iloc[0].copy()
    ret_merge["source"] = "ChimerAnno"
    ret_merge["type"] = args.new_type_name
    ret_merge["score"] = "."
    ret_merge["attributes"] = ""
    df.sort_values(["start", "end"], inplace=True, ascending=[True, False])
    df.drop_duplicates(inplace=True)
    df.reset_index(inplace=True, drop=True)
    for i in df.index.tolist()[:-1]:
        if is_reversed:
            bool_test = df.at[i, "source"] == "downstream" and df.at[i + 1, "source"] == "upstream"
        else:
            bool_test = df.at[i, "source"] == "upstream" and df.at[i + 1, "source"] == "downstream"
        if bool_test:
            ret_merge["start"] = df.at[i, "start"]
            ret_merge["end"] = df.at[i + 1, "end"]
            ret_merge["len"] = ret_merge["end"] - ret_merge["start"] + 1
            ret_df = ret_df.append(ret_merge, ignore_index=True)
    ret_df = ret_df[ret_df["len"] <= args.window_size]
    return ret_df


def _merge_perfect_overlaps_by_cluster(df: pd.DataFrame, args) -> pd.Series or pd.DataFrame:
    ret_merge = df.iloc[0].copy()
    ret_merge["source"] = "ChimerAnno"
    ret_merge["type"] = args.new_type_name
    ret_merge["score"] = "."
    ret_merge["attributes"] = ""
    is_reversed = True if "-" in df["strand"].tolist() else False
    upstream_start = df[df["source"] == "upstream"]["start"].min()
    upstream_end = df[df["source"] == "upstream"]["end"].max()
    downstream_start = df[df["source"] == "downstream"]["start"].min()
    downstream_end = df[df["source"] == "downstream"]["end"].max()
    if is_reversed and upstream_start >= downstream_start:
        ret_merge["start"] = downstream_start
        ret_merge["end"] = upstream_end
        return ret_merge
    elif not is_reversed and downstream_end >= upstream_end:
        ret_merge["start"] = upstream_start
        ret_merge["end"] = downstream_end
        return ret_merge
    else:
        return df


def _extend_downstream_ends(df: pd.DataFrame, args) -> pd.DataFrame:
    print("Adding offsets")
    up_df = df[df["source"] == "upstream"].copy()
    df.drop(up_df.index, inplace=True)
    up_df["len"] = up_df["end"] - up_df["start"] + 1
    up_df["offset"] = args.window_size - up_df["len"]
    f_up_df = up_df[up_df["strand"] == "+"].copy()
    r_up_df = up_df[up_df["strand"] == "-"].copy()
    del up_df
    f_up_df["original_start"], f_up_df["original_end"] = f_up_df["start"], f_up_df["end"]
    r_up_df["original_start"], r_up_df["original_end"] = r_up_df["start"], r_up_df["end"]
    f_up_df["end"] += f_up_df["offset"]
    r_up_df["start"] -= r_up_df["offset"]
    r_up_df.loc[r_up_df['start'] < 1, 'start'] = 1
    df = df.append(f_up_df, ignore_index=True)
    df = df.append(r_up_df, ignore_index=True)
    df.fillna(0, inplace=True)
    df[["len", "offset", "original_start", "original_end"]] = \
        df[["len", "offset", "original_start", "original_end"]].astype(int)
    return df


def _inject_prefix(in_path, prefix):
    dir = os.path.dirname(in_path)
    fname = os.path.basename(in_path)
    out_path = f"{dir}/{prefix}_{fname}" if dir != "" else f"{prefix}_{fname}"
    return out_path

main()