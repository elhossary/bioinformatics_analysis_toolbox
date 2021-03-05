import pandas as pd
from os import path
import argparse
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--main_gff", required=True, help="", type=str)
    parser.add_argument("--fetch_gffs", required=True, help="", type=str, nargs="+")
    parser.add_argument("--anno_type", required=True, help="", type=str)
    parser.add_argument("--stream", required=True, help="", type=str, choices=["up", "down"])
    parser.add_argument("--copy_attr", default=['all'], nargs="+")
    parser.add_argument("--copy_attr_like", default=None, help="", type=str, nargs="+")
    parser.add_argument("--out_gff", required=True, help="", type=str)
    args = parser.parse_args()
    args.copy_attr = [i.lower() for i in args.copy_attr]
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    main_gff_df = pd.read_csv(path.abspath(args.main_gff), names=col_names, sep="\t", comment="#")
    fetch_gffs_df = pd.DataFrame(columns=col_names)
    fetch_gffs_pathes = []
    for item in args.fetch_gffs:
        fetch_gffs_pathes.extend(glob.glob(item))

    for gff_path in fetch_gffs_pathes:
        x_gff_df = pd.read_csv(path.abspath(gff_path), names=col_names, sep="\t", comment="#")
        fetch_gffs_df = fetch_gffs_df.append(x_gff_df.copy())
        del x_gff_df

    fetch_gffs_df.drop_duplicates(inplace=True)
    fetch_gffs_df = fetch_gffs_df[fetch_gffs_df["type"] == args.anno_type].copy()
    for indx in main_gff_df.index:
        anno_start = main_gff_df.at[indx, 'start']
        anno_end = main_gff_df.at[indx, 'end']
        strand = main_gff_df.at[indx, 'strand']
        if strand == "+":
            pos = anno_end if args.stream == "down" else anno_start
            if args.stream == "down":
                x_df = fetch_gffs_df[(fetch_gffs_df['seqid'] == main_gff_df.at[indx, 'seqid']) &
                                     (fetch_gffs_df['strand'] == strand) &
                                     (fetch_gffs_df['start'] > pos)].copy()
                x_df.sort_values(["start"], inplace=True)
            else:
                x_df = fetch_gffs_df[(fetch_gffs_df['seqid'] == main_gff_df.at[indx, 'seqid']) &
                                     (fetch_gffs_df['strand'] == strand) &
                                     (fetch_gffs_df['end'] < pos)].copy()
                x_df.sort_values(["end"], inplace=True, ascending=False)
        elif strand == "-":
            pos = anno_end if args.stream == "up" else anno_start
            if args.stream == "up":
                x_df = fetch_gffs_df[(fetch_gffs_df['seqid'] == main_gff_df.at[indx, 'seqid']) &
                                     (fetch_gffs_df['strand'] == strand) &
                                     (fetch_gffs_df['start'] > pos)].copy()
                x_df.sort_values(["start"], inplace=True)
            else:
                x_df = fetch_gffs_df[(fetch_gffs_df['seqid'] == main_gff_df.at[indx, 'seqid']) &
                                     (fetch_gffs_df['strand'] == strand) &
                                     (fetch_gffs_df['end'] < pos)].copy()
                x_df.sort_values(["end"], inplace=True, ascending=False)
        else:
            x_df = None
            exit(1)

        if x_df.empty:
            continue
        x_first_row = x_df.iloc[0, :]
        x_attr = {}

        dist = x_first_row["start"] - anno_end if strand == "+" else anno_start - x_first_row["end"]
        x_attr[f"{args.stream}stream_{x_first_row['type']}_distance"] = str(dist)
        x_first_row_attr = parse_attributes(x_first_row['attributes'])
        for k, v in x_first_row_attr.items():
            if not args.copy_attr == ["all"] and k.lower() not in args.copy_attr:
                continue
            x_attr[f"{args.stream}stream_{x_first_row['type']}_{k}"] = v
        for k, v in x_first_row_attr.items():
            for attr_like in args.copy_attr_like:
                if attr_like not in k:
                    continue
                x_attr[f"{args.stream}stream_{x_first_row['type']}_{k}"] = v
        main_gff_df.at[indx, "attributes"] += f";{gen_attr(x_attr)}"
        main_gff_df.at[indx, "attributes"] = main_gff_df.at[indx, "attributes"].replace(";;", ";").strip(";")
    print("Writing GFF file...")
    main_gff_df.to_csv(path.abspath(args.out_gff), sep="\t", header=False, index=False)


def gen_attr(dict_in):
    out_str = ""
    for k, v in dict_in.items():
        out_str += f"{k}={v};"
    return out_str[:-1]


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}


main()
exit(0)
