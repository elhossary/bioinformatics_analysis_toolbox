import pandas as pd
import os
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tab", required=True, help="", type=str)
    parser.add_argument("--link", required=True, help="", type=str)
    parser.add_argument("--tracks", required=True, help="", type=str)
    args = parser.parse_args()
    tab_df = pd.read_excel(os.path.abspath(args.tab))
    tab_df["Browser_link"] = ""
    for i in tab_df.index:
        tab_df.at[i, "Browser_link"] = generate_link(args.link, args.tracks, tab_df.iloc[i, :])

    tab_df.to_excel(os.path.abspath(args.tab), header=True, index=False)


def generate_link(url, tracks, row):
    output_str = f"{url}&loc={row['seq_id']}%3A" \
                  f"{str(int(row['start']) - 30) if int(row['start']) - 30 > 0 else 0}" \
                  f"..{str(int(row['end']) + 30)}" \
                  f"&tracks={tracks}"
    return output_str


main()
