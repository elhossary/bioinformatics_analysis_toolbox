import argparse
import pandas as pd
from os import path
parser = argparse.ArgumentParser()
parser.add_argument("--gff_in", required=True, help="", type=str)
parser.add_argument("--poly_t_in", required=True, help="", type=str)
parser.add_argument("--gff_out", required=True, help="", type=str)
args = parser.parse_args()

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
poly_t_df = pd.read_csv(path.abspath(args.poly_t_in), names=col_names, sep="\t", comment="#")
for index, row in gff_df.iterrows():
    if row['strand'] == "+":
        tmp = poly_t_df[(poly_t_df['seqid'] == row['seqid']) & (poly_t_df['strand'] == "+") &
                        (((poly_t_df['start'] <= row['end']) &
                         (row['end'] <= poly_t_df['end'])))]
        if tmp.empty:
            tmp = poly_t_df[(poly_t_df['seqid'] == row['seqid']) & (poly_t_df['strand'] == "+")
                            & (poly_t_df['start'].between(row['start'], row['end'], inclusive=True))
                            & (poly_t_df['end'].between(row['start'], row['end'], inclusive=True))]

        if tmp.empty:
            print(tmp)
            print(row['end'])
        #print(tmp.to_string())