import pandas as pd
from glob import glob
import argparse
import sys
import os


def main():
    # Params
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_path", required=True, help="ANNOgesic project(s) folder path", type=str)

    args = parser.parse_args()
    # ---------------------------
    folders_list = glob(args.project_path)
    if len(folders_list) <= 0:
        print("Error")
        sys.exit(os.EX_USAGE)

main()
