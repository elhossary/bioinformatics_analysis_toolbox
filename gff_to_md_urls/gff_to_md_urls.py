import pandas as pd
import argparse
import glob
import os
import datetime
import glob

def main():
    # Param
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="GFF(s) file input", type=str)
    parser.add_argument("--url", required=True, help="URL to the JBrowse folder", type=str)
    parser.add_argument("--tracks", required=False, help="", type=str, default="")
    parser.add_argument("--extra_column", required=False, help="", type=str, nargs='+')
    parser.add_argument("--checkin_files", required=False, help="", type=str, nargs='+')
    args = parser.parse_args()
    gff_files = glob.glob(args.gff_in)
    extra_columns = ""
    extra_columns_sep = ""
    checks_list = []

    if args.extra_column is not None:
        for col in args.extra_column:
            extra_columns += f"|**{col}**"
            extra_columns_sep += f"|:-----:"
        for item in args.checkin_files:
            checks_list.append(open(os.path.abspath(item), "r").read())

    column_names = ["accession", "source", "type", "start", "end", "dot1", "strand", "dot2", "attributes"]
    header = f"**No.**|**Name**|**Genomic Location**{extra_columns}\n" \
             f":-----:|:-----:|:-----:{extra_columns_sep}\n"
    output_str = ""
    counter = 0
    for file in gff_files:
        gff_df = pd.read_csv(os.path.abspath(file), names=column_names, sep="\t", comment="#")
        for index, row in gff_df.iterrows():
            counter += 1
            output_str += f"{counter}|[{get_label_name(parse_attributes(row['attributes']))}]" \
                          f"({args.url}&loc={row['accession']}%3A" \
                          f"{str(int(row['start']) - 30) if int(row['start']) - 30 > 0 else 0}" \
                          f"..{str(int(row['end']) + 30)}" \
                          f"&highlight={row['accession']}%3A{row['start']}..{row['end']}" \
                          f"&tracks={args.tracks})|" \
                          f"{row['accession']} .. {row['start']} .. {row['end']} .. {row['strand']}" \
                          f"{check(checks_list, row['attributes'])}" \
                          f"\n"
        # {extra_columns_sep.replace(':-----:', '')}
        counter = 0
        output_str = header + output_str
        output_path = os.path.abspath(os.path.join(file, os.pardir))
        output_basename = os.path.basename(file).replace('.gff', '.md').replace('.GFF', '.md')
        out_file = open(f"{output_path}/{datetime.datetime.today().strftime('%Y-%m-%d')}_{output_basename}", "w")
        out_file.write(output_str)
        out_file.close()


def get_label_name(dict_in):
    if 'name' in dict_in.keys():
        return dict_in['name']
    elif 'label' in dict_in.keys():
        return dict_in['label']
    elif 'id' in dict_in.keys():
        return dict_in['id']
    else:
        return 'Click here'

def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

def check(checks_list, attr):
    checks_str = ""
    for check in checks_list:
        if get_label_name(parse_attributes(attr)) in check:
            checks_str += f"|Yes"
        else:
            checks_str += f"|No"
    return checks_str
main()