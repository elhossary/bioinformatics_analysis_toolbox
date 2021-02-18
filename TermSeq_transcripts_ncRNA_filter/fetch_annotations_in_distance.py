import pandas as pd
from os import path
import argparse
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--main_gff", required=True, help="", type=str)
    parser.add_argument("--fetch_gffs", required=True, help="", type=str, nargs="+")
    parser.add_argument("--stream", required=True, help="", type=str, choices=["up", "down"])
    parser.add_argument("--distance", default=50, required=False, help="", type=int)
    parser.add_argument("--copy_attr", default=['all'], nargs="+")
    parser.add_argument("--copy_columns", default=None,
                        choices=["seqid", "source", "type", "start", "end", "score", "strand", "phase"], nargs="+")
    parser.add_argument("--allow_overlap", default=False, action='store_true')
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

    for indx in main_gff_df.index:

        anno_start = main_gff_df.at[indx, 'start']
        anno_end = main_gff_df.at[indx, 'end']
        strand = main_gff_df.at[indx, 'strand']
        anno_range = range(anno_start, anno_end + 1, 1)
        if strand == "+":
            if args.allow_overlap:
                search_start = anno_start if args.stream == "down" else anno_start - 1 - args.distance
                search_end = anno_end + 1 + args.distance if args.stream == "down" else anno_end
            else:
                search_start = anno_end + 1 if args.stream == "down" else anno_start - 1 - args.distance
                search_end = anno_end + 1 + args.distance if args.stream == "down" else anno_start - 1
            x_df = fetch_gffs_df[(fetch_gffs_df['seqid'] == main_gff_df.at[indx, 'seqid']) &
                                 (fetch_gffs_df['strand'] == strand) &
                                 (fetch_gffs_df['start'].isin(range(search_start, search_end + 2, 1)))].copy()
        elif strand == "-":
            if args.allow_overlap:
                search_start = anno_start if args.stream == "up" else anno_start - 1 - args.distance
                search_end = anno_end + 1 + args.distance if args.stream == "up" else anno_end
            else:
                search_start = anno_end + 1 if args.stream == "up" else anno_start - 1 - args.distance
                search_end = anno_end + 1 + args.distance if args.stream == "up" else anno_start - 1
            x_df = fetch_gffs_df[(fetch_gffs_df['seqid'] == main_gff_df.at[indx, 'seqid']) &
                                 (fetch_gffs_df['strand'] == strand) &
                                 (fetch_gffs_df['end'].isin(range(search_start, search_end + 2, 1)))].copy()
        else:
            x_df = None
            exit(1)

        if x_df.empty:
            continue
        x_df.sort_values(['type', 'start', 'end'], inplace=True)
        x_attr = {}
        for x_indx in x_df.index:
            dist = x_df.at[x_indx, "start"] - anno_end - 1 if strand == "+" else anno_start - x_df.at[x_indx, "end"]- 1
            if f"{args.stream}stream_{x_df.at[x_indx, 'type']}_distance" not in x_attr.keys():
                x_attr[f"{args.stream}stream_{x_df.at[x_indx, 'type']}_distance"] = str(dist)
            else:
                x_attr[f"{args.stream}stream_{x_df.at[x_indx, 'type']}_distance"] += f"|{str(dist)}"

            for k, v in parse_attributes(x_df.at[x_indx, 'attributes']).items():
                if not args.copy_attr == ["all"] and k.lower() not in args.copy_attr:
                    continue
                if f"{args.stream}stream_{x_df.at[x_indx, 'type']}_{k}" not in x_attr.keys():
                    x_attr[f"{args.stream}stream_{x_df.at[x_indx, 'type']}_{k}"] = v
                else:
                    x_attr[f"{args.stream}stream_{x_df.at[x_indx, 'type']}_{k}"] += f"|{v}"

            if args.copy_columns is not None:
                for col in args.copy_columns:
                    if f"{args.stream}stream_{x_df.at[x_indx, 'type']}_{col}" not in x_attr.keys():
                        x_attr[f"{args.stream}stream_{x_df.at[x_indx, 'type']}_{col}"] = x_df.at[x_indx, col]
                    else:
                        x_attr[f"{args.stream}stream_{x_df.at[x_indx, 'type']}_{col}"] += f"|{x_df.at[x_indx, col]}"

            is_overlapping = False
            if args.allow_overlap:
                if x_df.at[x_indx, 'start'] in anno_range or x_df.at[x_indx, 'start'] in anno_range:
                    is_overlapping = True
                if f"{args.stream}stream_{x_df.at[x_indx, 'type']}_is_overlapping" not in x_attr.keys():
                    x_attr[f"{args.stream}stream_{x_df.at[x_indx, 'type']}_is_overlapping"] = \
                        'YES' if is_overlapping else 'NO'
                else:
                    x_attr[f"{args.stream}stream_{x_df.at[x_indx, 'type']}_is_overlapping"] += \
                        f"|{'YES' if is_overlapping else 'NO'}"
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
