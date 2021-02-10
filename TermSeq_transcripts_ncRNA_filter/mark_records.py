import pandas as pd
import os
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tab1", required=True, help="", type=str)
    parser.add_argument("--tab2", required=True, help="", type=str)
    args = parser.parse_args()

    tab1_df = pd.read_excel(os.path.abspath(args.tab1))
    tab2_df = pd.read_excel(os.path.abspath(args.tab2))
    tab1_df.fillna("", inplace=True)
    tab2_df.fillna("", inplace=True)
    col_name = "downstream CDS e coli k12 mg1655 ortholog protein id"
    f_col_name = "protein id found in k. pneumonia table"
    f_id_col_name = "protein id name in k. pneumonia table"
    tab1_df = find_commons(tab1_df, tab2_df, col_name, f_col_name, f_id_col_name)
    col_name = "downstream CDS k pneumoniae kpnih1 ortholog protein id"
    f_col_name = "protein id found in e. coli table"
    f_id_col_name = "protein id name in e. coli table"
    tab2_df = find_commons(tab2_df, tab1_df, col_name, f_col_name, f_id_col_name)

    tab1_path = os.path.dirname(args.tab1)
    tab1_file_name = os.path.basename(args.tab1)
    tab1_save_path = f"{tab1_path}/updated_{tab1_file_name}"
    tab1_df.to_excel(os.path.abspath(tab1_save_path), index=False)

    tab2_path = os.path.dirname(args.tab2)
    tab2_file_name = os.path.basename(args.tab2)
    tab2_save_path = f"{tab2_path}/updated_{tab2_file_name}"
    tab2_df.to_excel(os.path.abspath(tab2_save_path), index=False)


def find_commons(tab1_df, tab2_df, col_name, f_col_name, f_id_col_name):
    tab1_df[f_col_name] = "NO"
    tab1_df[f_id_col_name] = ""
    for x in tab1_df.index:
        find_word = tab1_df.at[x, "downstream CDS name"]
        if find_word != "":
            f_df = tab2_df[tab2_df[col_name].str.contains(find_word)]
            if f_df.empty:
                continue
            tab1_df.at[x, f_col_name] = "YES"
            tab1_df.at[x, f_id_col_name] = ','.join(f_df["downstream CDS name"].values.tolist())


    return tab1_df


main()