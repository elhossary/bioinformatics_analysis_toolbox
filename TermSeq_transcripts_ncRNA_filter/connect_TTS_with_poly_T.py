import argparse
import pandas as pd
import os
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--tts_in", required=True, help="", type=str)
parser.add_argument("--poly_t_in", required=True, help="", type=str)
parser.add_argument("--range", required=True, help="", type=int)
parser.add_argument("--gff_out", required=True, help="", type=str)
args = parser.parse_args()

col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
tts_df = pd.read_csv(os.path.abspath(args.tts_in), names=col_names, sep="\t", comment="#")
poly_t_df = pd.read_csv(os.path.abspath(args.poly_t_in), names=col_names, sep="\t", comment="#")
out_list = []
for idx in tts_df.index:
    if tts_df.at[idx, "strand"] == "+":
        start = tts_df.at[idx, "start"] - args.range
        if start < 1:
            start = 1
        end = tts_df.at[idx, "start"] + args.range
        # select all possiblites
        tmp_df = poly_t_df[(poly_t_df["seqid"] == tts_df.at[idx, "seqid"]) &
                           (poly_t_df["strand"] == tts_df.at[idx, "strand"]) &
                           ((poly_t_df["start"].between(start, end)) | (poly_t_df["end"].between(start, end)))]
        # choose the nearset
        poly_start = 0
        poly_end = 0
        if not tmp_df.empty:
            if tmp_df.shape[0] > 1:

                print(tmp_df.to_string())
                poly_start = 0
                poly_end = 0
            else:
                poly_start = tmp_df["start"].iloc[0]
                poly_end = tmp_df["end"].iloc[0]
        if poly_start <= tts_df.at[idx, "start"] <= poly_end:
            out_list.append([poly_start, poly_end, "+"])
        elif poly_start > tts_df.at[idx, "start"]:
            out_list.append([tts_df.at[idx, "start"], poly_end, "+"])
        elif poly_end < tts_df.at[idx, "start"]:
            out_list.append([poly_start, tts_df.at[idx, "start"], "+"])
        else:
            print(f"Fatal error, {tts_df.at[idx, 'start']} - {poly_start} - {poly_end}")
    elif tts_df.at[idx, "strand"] == "-":
        start = tts_df.at[idx, "end"] + args.range
        if start < 1:
            start = 1
        end = tts_df.at[idx, "end"] - args.range
        # select all possiblites
        tmp_df = poly_t_df[(poly_t_df["start"].between(start, end)) | (poly_t_df["end"].between(start, end))]
        print(tmp_df)
    else:
        print("Fatal error")
        exit(1)