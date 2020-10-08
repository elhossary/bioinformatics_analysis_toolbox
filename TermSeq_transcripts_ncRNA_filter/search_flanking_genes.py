import argparse
import pandas as pd
from os import path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--ref_gff_in", required=True, help="", type=str)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    parser.add_argument("--force_strandedness", default=False, help="", action='store_true')
    parser.add_argument("--combine_flanks_info", default=False, help="", action='store_true')
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    ref_gff_df = pd.read_csv(path.abspath(args.ref_gff_in), names=col_names, sep="\t", comment="#")
    ref_gff_df = ref_gff_df[ref_gff_df["type"] == "gene"]

    df_dict = {}
    for seqid in gff_df["seqid"].unique():
        df_dict[seqid] = ref_gff_df[ref_gff_df["seqid"] == seqid]
        if df_dict[seqid].empty:
            df_dict[seqid] = None
            print(f"Cannot find {seqid} in the reference annotations")
    gff_df = gff_df.apply(lambda x: append_flanking_genes_to_attributes(x, df_dict[x.seqid], args.force_strandedness,
                                                                        args.combine_flanks_info), axis=1)
    print("Writing GFF file")
    gff_df.to_csv(path.abspath(args.gff_out), sep="\t", index=False, header=False)


def append_flanking_genes_to_attributes(gff_row, ref_df, strandedness, combine):
    if ref_df is None:
        return gff_row
    downstream_distance = 0
    downstream_gene = ""
    upstream_distance = 0
    upstream_gene = ""
    downstream_gene_strand = ""
    upstream_gene_strand = ""
    overlap_rows = None

    if gff_row.strand == "+":
        f_ref_df = ref_df
        overlap_rows = f_ref_df[((f_ref_df["start"].isin(range(gff_row.start + 1, gff_row.end, 1))) &
                                 (f_ref_df["strand"] == "+")) |
                                ((f_ref_df["end"].isin(range(gff_row.start + 1, gff_row.end, 1))) &
                                 (f_ref_df["strand"] == "+"))]
        if strandedness:
            f_ref_df = ref_df[ref_df["strand"] == "+"]
        try:
            up_row = f_ref_df[f_ref_df["end"] <= gff_row.start].sort_values(["end"], ascending=False).iloc[0]
            upstream_gene = parse_attributes(up_row["attributes"])["name"]
            upstream_distance = gff_row.start - up_row["end"]
            if up_row.strand == "+":
                upstream_gene_strand = "sense"
            else:
                upstream_gene_strand = "antisense"
        except:
            upstream_gene = "NONE"
            upstream_distance = "NONE"
            upstream_gene_strand = "NONE"
        try:
            down_row = f_ref_df[f_ref_df["start"] >= gff_row.end].sort_values(["start"]).iloc[0]
            downstream_gene = parse_attributes(down_row["attributes"])["name"]
            downstream_distance = down_row["start"] - gff_row.end
            if down_row.strand == "+":
                downstream_gene_strand = "sense"
            else:
                downstream_gene_strand = "antisense"
        except:
            downstream_gene = "NONE"
            downstream_distance = "NONE"
            downstream_gene_strand = "NONE"

    elif gff_row.strand == "-":
        r_ref_df = ref_df
        overlap_rows = r_ref_df[((r_ref_df["start"].isin(range(gff_row.start, gff_row.end + 1, 1))) &
                                 (r_ref_df["strand"] == "-")) |
                                ((r_ref_df["end"].isin(range(gff_row.start, gff_row.end + 1, 1))) &
                                 (r_ref_df["strand"] == "-"))]
        if strandedness:
            r_ref_df = ref_df[ref_df["strand"] == "-"]
        try:
            down_row = r_ref_df[r_ref_df["end"] <= gff_row.start].sort_values(["end"], ascending=False).iloc[0]
            downstream_gene = parse_attributes(down_row["attributes"])["name"]
            downstream_distance = gff_row.start - down_row["end"]
            if down_row.strand == "-":
                downstream_gene_strand = "sense"
            else:
                downstream_gene_strand = "antisense"
        except:
            downstream_gene = "NONE"
            downstream_distance = "NONE"
            downstream_gene_strand = "NONE"
        try:
            up_row = r_ref_df[r_ref_df["start"] >= gff_row.end].sort_values(["start"]).iloc[0]
            upstream_gene = parse_attributes(up_row["attributes"])["name"]
            upstream_distance = up_row["start"] - gff_row.end
            if up_row.strand == "-":
                upstream_gene_strand = "sense"
            else:
                upstream_gene_strand = "antisense"
        except:
            upstream_gene = "NONE"
            upstream_distance = "NONE"
            upstream_gene_strand = "NONE"

    else:
        print("Fatal Error")
        exit(1)

    if not overlap_rows.empty and overlap_rows is not None:
        overlapping_genes = '|'.join([parse_attributes(i)["name"] for i in overlap_rows["attributes"].values.tolist()])
    else:
        overlapping_genes = "NONE"
    if combine:
        gff_row.attributes += f";up_flank={upstream_gene}|{upstream_gene_strand}|{upstream_distance}" \
                              f";down_flank={downstream_gene}|{downstream_gene_strand}|{downstream_distance}" \
                              f";overlapping_genes={overlapping_genes}"
    else:
        gff_row.attributes += f";up_flank={upstream_gene}" \
                              f";up_flank_dist={upstream_distance}" \
                              f";up_flank_strand={upstream_gene_strand}" \
                              f";down_flank={downstream_gene}" \
                              f";down_flank_dist={downstream_distance}" \
                              f";down_flank_strand={downstream_gene_strand}" \
                              f";overlapping_genes={overlapping_genes}"
    return gff_row

def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()
