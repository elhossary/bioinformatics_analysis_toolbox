import argparse
import pandas as pd
from os import path


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--ref_gff_in", required=True, help="", type=str)
parser.add_argument("--gff_out", required=True, help="", type=str)
args = parser.parse_args()
col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
ref_df = pd.read_csv(path.abspath(args.ref_gff_in), names=col_names, sep="\t", comment="#")
str_out = ""
for index, gff_row in gff_df.iterrows():
    tmp_str = ""
    downstream_distance = 0
    downstream_gene = ""
    upstream_distance = 0
    upstream_gene = ""
    downstream_gene_strand = ""
    upstream_gene_strand = ""
    overlapping_genes = ""
    tmp = None
    f_tmp = None
    r_tmp = None

    if gff_row["strand"] == "+":
        tmp = ref_df[(ref_df["seqid"] == gff_row["seqid"]) & (ref_df["type"] == "gene") &
                     (ref_df["start"] >= gff_row["end"])].sort_values(["start"]).head(1)
        if not tmp.empty:
            downstream_gene = parse_attributes(tmp.iloc[0]['attributes'])["name"]
            downstream_distance = tmp.iloc[0]['start'] - gff_row["end"]
            if tmp.iloc[0]['strand'] == "+":
                downstream_gene_strand = "same"
            else:
                downstream_gene_strand = "opposite"

        tmp = ref_df[(ref_df["seqid"] == gff_row["seqid"]) & (ref_df["type"] == "gene") &
                     (ref_df["end"] <= gff_row["start"])].sort_values(["end"], ascending=False).head(1)
        if not tmp.empty:
            upstream_gene = parse_attributes(tmp.iloc[0]['attributes'])["name"]
            upstream_distance = gff_row["start"] - tmp.iloc[0]['end']
            if tmp.iloc[0]['strand'] == "+":
                upstream_gene_strand = "same"
            else:
                upstream_gene_strand = "opposite"
        tmp = ref_df[(ref_df["seqid"] == gff_row["seqid"]) & (ref_df["strand"] == "+") & (ref_df["type"] == "gene") &
                     (((ref_df["start"] <= gff_row["start"]) & (gff_row["start"] <= ref_df["end"])) |
                      ((ref_df["start"] <= gff_row["end"]) & (gff_row["end"] <= ref_df["end"])))]
        if not tmp.empty:
            for name in [parse_attributes(x)["name"] for x in tmp["attributes"].values.tolist()]:
                if overlapping_genes == "":
                    overlapping_genes += name
                else:
                    overlapping_genes += f",{name}"
        else:
            overlapping_genes = "#"
    elif gff_row["strand"] == "-":
        
        tmp = ref_df[(ref_df["seqid"] == gff_row["seqid"]) & (ref_df["type"] == "gene") &
                     (ref_df["end"] <= gff_row["start"])].sort_values(["end"], ascending=False).head(1)
        if not tmp.empty:
            downstream_gene = parse_attributes(tmp.iloc[0]['attributes'])["name"]
            downstream_distance = gff_row["start"] - tmp.iloc[0]['end']
            if tmp.iloc[0]['strand'] == "-":
                downstream_gene_strand = "same"
            else:
                downstream_gene_strand = "opposite"
        tmp = ref_df[(ref_df["seqid"] == gff_row["seqid"]) & (ref_df["type"] == "gene") &
                     (ref_df["start"] >= gff_row["end"])].sort_values(["start"]).head(1)

        if not tmp.empty:
            upstream_gene = parse_attributes(tmp.iloc[0]['attributes'])["name"]
            upstream_distance = tmp.iloc[0]['start'] - gff_row["end"]
            if tmp.iloc[0]['strand'] == "-":
                upstream_gene_strand = "same"
            else:
                upstream_gene_strand = "opposite"
        tmp = ref_df[(ref_df["seqid"] == gff_row["seqid"]) & (ref_df["strand"] == "-") & (ref_df["type"] == "gene") &
                     (((ref_df["start"] <= gff_row["start"]) & (gff_row["start"] <= ref_df["end"])) |
                      ((ref_df["start"] <= gff_row["end"]) & (gff_row["end"] <= ref_df["end"])))]
        if not tmp.empty:
            for name in [parse_attributes(x)["name"] for x in tmp["attributes"].values.tolist()]:
                if overlapping_genes == "":
                    overlapping_genes += name
                else:
                    overlapping_genes += f",{name}"
        else:
            overlapping_genes = "#"
    else:
        print("Fatal error")
        exit()

    tmp_str = f";up_gene={upstream_gene};up_gene_dist={upstream_distance};up_gene_strand={upstream_gene_strand}" \
              f";down_gene={downstream_gene};down_gene_dist={downstream_distance};down_gene_strand={downstream_gene_strand}" \
              f";overlapping_genes={overlapping_genes}"
    str_out += \
        f"{gff_row['seqid']}\t" + \
        f"{gff_row['source']}\t" + \
        f"{gff_row['type']}\t" + \
        f"{int(gff_row['start'])}\t" + \
        f"{int(gff_row['end'])}\t" + \
        f"{gff_row['score']}\t" + \
        f"{gff_row['strand']}\t" + \
        f"{gff_row['phase']}\t" + \
        f"{gff_row['attributes']}{tmp_str}" + \
        "\n"
print("Writing GFF file...")
outfile = open(path.abspath(args.gff_out), "w")
outfile.write(f"{str_out}")
outfile.close()
