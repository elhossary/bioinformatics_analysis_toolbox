import argparse
import pandas as pd
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True, help="", type=str)
    parser.add_argument("--bases", required=True, help="", type=int)
    parser.add_argument("--prime_end", required=True, help="", type=str, choices=["both", "5", "3"], default="both")
    args = parser.parse_args()

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_file = os.path.abspath(args.gff)
    df = pd.read_csv(gff_file, sep="\t", names=col_names, comment="#")

    if args.prime_end == "both":
        df["start"] += args.bases
        df["end"] -= args.bases
    else:
        f_df = df[df["strand"] == "+"].copy()
        r_df = df[df["strand"] == "-"].copy()
        if args.prime_end == "5":
            f_df["start"] += args.bases
            r_df["end"] -= args.bases
        elif args.prime_end == "3":
            f_df["end"] -= args.bases
            r_df["start"] += args.bases
        else:
            exit(1)
        df = f_df.append(r_df)
    df["len"] = df["end"] - df["start"]
    df = df[df["len"] > 0]
    df.drop(["len"], axis=1, inplace=True)
    df.sort_values(["seqid", "start", "end"], inplace=True)
    df.to_csv(gff_file, sep="\t", header=False, index=False)
    print("Annotations trimmed")


if __name__ == '__main__':
    main()