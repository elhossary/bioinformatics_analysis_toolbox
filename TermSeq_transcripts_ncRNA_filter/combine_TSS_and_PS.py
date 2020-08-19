import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tss_in", required=True, help="", type=str)
    parser.add_argument("--ps_in", required=True, help="", type=str)
    parser.add_argument("--tss_ps_out", required=True, help="", type=str)
    parser.add_argument("--split", default=False, help="", action='store_true')
    args = parser.parse_args()
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    tss_df = pd.read_csv(os.path.abspath(args.tss_in), names=col_names, sep="\t", comment="#")
    ps_df = pd.read_csv(os.path.abspath(args.ps_in), names=col_names, sep="\t", comment="#")
    out_df = pd.DataFrame(columns=col_names)
    out_df = out_df.append(tss_df, ignore_index=True)
    out_df = out_df.append(ps_df, ignore_index=True)
    out_df.sort_values(["seqid", "start", "end", "strand", "type"], inplace=True)
    out_df.reset_index(inplace=True, drop=True)
    out_df.drop_duplicates(subset=["seqid", "start", "end", "strand"], keep="first", inplace=True)
    tss_ps_name = str(os.path.abspath(f"{args.tss_ps_out}"))
    out_df.to_csv(tss_ps_name, sep="\t", header=False, index=False)
    if args.split:
        print(tss_ps_name.replace('_PS', ''))
        tss_name = tss_ps_name.replace('_PS', '')
        ps_name = tss_ps_name.replace('_TSS', '')
        tss_out_df = out_df[out_df["type"] == "TSS"]
        ps_out_df = out_df[out_df["type"] == "processing_site"]
        tss_out_df.to_csv(tss_name, sep="\t", header=False, index=False)
        ps_out_df.to_csv(ps_name, sep="\t", header=False, index=False)

main()