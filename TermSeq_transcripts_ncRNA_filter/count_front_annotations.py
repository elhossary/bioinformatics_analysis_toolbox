import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in_a", required=True, help="", type=str)
    parser.add_argument("--gff_in_b", required=True, help="", type=str)
    parser.add_argument("--range", default=50, help="", type=int)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df_a = pd.read_csv(args.gff_in_a, names=col_names, sep="\t", comment="#")
    gff_df_b = pd.read_csv(args.gff_in_b, names=col_names, sep="\t", comment="#")
    for indx_a in gff_df_a.index:
        if gff_df_a.at[indx_a, "strand"] == "+":
            x = gff_df_b[(gff_df_b["seqid"] == gff_df_a.at[indx_a, "seqid"]) and
                          (gff_df_b["strand"] == "+") and
                          (gff_df_b["end"].between(gff_df_a.at[indx_a, "end"], gff_df_a.at[indx_a, "end"] + args.range))]
        elif gff_df_a.at[indx_a, "strand"] == "-":
            x = gff_df_b[(gff_df_b["seqid"] == gff_df_a.at[indx_a, "seqid"]) and
                         (gff_df_b["strand"] == "-") and
                         (gff_df_b["start"].between(gff_df_a.at[indx_a, "start"] - args.range, gff_df_a.at[indx_a, "start"]))]
        else:
            exit(1)
        print(x)
main()