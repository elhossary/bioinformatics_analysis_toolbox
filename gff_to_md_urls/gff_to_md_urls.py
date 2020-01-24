import pandas as pd
import argparse
import glob
import os
import datetime


def main():
    # Param
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="GFF(s) file input", type=str)
    parser.add_argument("--url", required=True, help="URL to the JBrowse folder", type=str)
    args = parser.parse_args()
    gff_files = glob.glob(args.gff_in)

    column_names = ["accession", "source", "type", "start", "end", "dot1", "strand", "dot2", "attributes"]

    output_str = ""
    for file in gff_files:
        gff_df = pd.read_csv(os.path.abspath(file), names=column_names, sep="\t", comment="#")
        for index, row in gff_df.iterrows():
            output_str += f"- [{row['accession']} {row['start']}:{row['end']} {row['strand']}]"\
                           f"({args.url}&loc={row['accession']}%3A"\
                           f"{str(int(row['start']) - 20) if int(row['start']) - 20 > 0 else 0}"\
                           f"..{str(int(row['end']) + 20)}"\
                           f"&highlight={row['accession']}%3A{row['start']}..{row['end']})\n"
        output_path = os.path.abspath(os.path.join(file, os.pardir))
        output_basename = os.path.basename(file).replace('.gff', '.md').replace('.GFF', '.md')
        out_file = open(f"{output_path}/{datetime.datetime.today().strftime('%Y-%m-%d')}_{output_basename}", "w")
        out_file.write(output_str)
        out_file.close()


main()

