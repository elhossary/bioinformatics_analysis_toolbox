from wiggle_parser.wiggle_parser import WiggleParser as wp
import argparse
from numpy import genfromtxt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--f_wig_in", required=True, help="", type=str)
    parser.add_argument("--r_wig_in", required=True, help="", type=str)
    parser.add_argument("--tss_gff_in", required=True, help="", type=str)
    parser.add_argument("--min_coverage", required=True, help="", type=float)
    args = parser.parse_args()
    accession_list = ['NC_002505.1', 'NC_002506.1']
    f_wig_parsed = wp(args.f_wig_in).parse()
    r_wig_parsed = wp(args.r_wig_in).parse()
    counter = 0
    for accession in accession_list:
        f_wig_df_sliced = f_wig_parsed[accession][f_wig_parsed[accession][1] >= args.min_coverage]
        r_wig_df_sliced = r_wig_parsed[accession][r_wig_parsed[accession][1] <= args.min_coverage * -1]
        tss_arr = build_arr_form_gff(args.tss_gff_in)
        for tss_index, tss_row in enumerate(tss_arr):
            if tss_row[0] == accession:
                if tss_row[6] == "+":
                    if get_score_of_wig_loc(f_wig_df_sliced, tss_row[3]):
                        counter += 1
                elif tss_row[6] == "-":
                    if get_score_of_wig_loc(r_wig_df_sliced, tss_row[3]):
                        counter += 1
                else:
                    print("ERROR")
                    exit()
    print(counter)


def get_score_of_wig_loc(wig_df, pos):
    x = wig_df[wig_df[0].isin(range(pos - 1, pos + 1))]
    if x.empty:
        return False
    else:
        return True


def build_arr_form_gff(path):
    data_arr = genfromtxt(path, delimiter="\t", comments="#", dtype=None, encoding=None)
    return data_arr

main()

