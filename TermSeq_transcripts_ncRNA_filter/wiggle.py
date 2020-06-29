import pandas as pd
import os
import numpy as np


class Wiggle:
    def __init__(self, file_path):
        self.file_path = file_path
        self.ret_df = None

    def parse(self):
        tmp_list = []
        current_wiggle_meta = {}
        column_names = ["track_type", "track_name", "variableStep_chrom", "variableStep_span", "location", "score"]
        with open(self.file_path, "r") as raw_file:
            print(f"==> Loaded file: {os.path.basename(self.file_path)}")
            for line in raw_file.readlines():
                if line[0].isnumeric():
                    tmp_list.append([current_wiggle_meta["track_type"],
                                     current_wiggle_meta["track_name"],
                                     current_wiggle_meta["variableStep_chrom"],
                                     current_wiggle_meta["variableStep_span"],
                                     int(line.split(" ")[0]),
                                     float(line.split(" ")[1].replace("\n", ""))])
                else:
                    current_wiggle_meta = self.parse_wiggle_header(line, current_wiggle_meta)
                    if len(current_wiggle_meta.keys()) == 2:
                        print(f"===> Parsing condition: {current_wiggle_meta['track_name']}")
                    elif len(current_wiggle_meta.keys()) == 4:
                        print(f"====> For sequence ID: {current_wiggle_meta['variableStep_chrom']}")
                    else:
                        exit(1)
        self.ret_df = pd.DataFrame(tmp_list, columns=column_names)
        return self.ret_df

    @staticmethod
    def parse_wiggle_header(line, current_wiggle_meta):
        if "type=" in line:
            current_wiggle_meta["track_type"] = line.split('type=')[-1].split(' ')[0].replace('\"', '')
        if "name=" in line:
            current_wiggle_meta["track_name"] = line.split('name=')[-1].replace('\n', '').replace('\"', '')
        if "chrom=" in line:
            current_wiggle_meta["variableStep_chrom"] = line.split('chrom=')[-1].split(' ')[0].replace('\"', '')
        if "span=" in line:
            current_wiggle_meta["variableStep_span"] = line.split('span=')[-1].replace('\n', '').replace('\"', '')
        return current_wiggle_meta

    def get_percentile(self, nth_percentile):
        if self.ret_df is None:
            print("Error: Wiggle isn't parsed")
            exit()
        return np.percentile(self.ret_df["score"].values.tolist(), nth_percentile)

    def normalize(self, factor):
        if self.ret_df is None:
            print("Error: Wiggle isn't parsed")
            exit()
        normalized_ret_df = self.ret_df
        normalized_ret_df.score = normalized_ret_df.score / factor
        return normalized_ret_df
