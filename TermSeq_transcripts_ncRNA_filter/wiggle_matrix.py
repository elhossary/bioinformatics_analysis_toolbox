from wiggle import Wiggle
import pandas as pd
import os
import glob
from Bio import SeqIO

class WiggleMatrix:
    def __init__(self, fasta_file, wiggle_files):
        self.fasta_file = fasta_file
        self.wiggle_files = wiggle_files
        self.wiggle_matrix_df = None

    def build_matrix(self):
        parsed_fasta = SeqIO.parse(self.fasta_file, "fasta")
        parsed_wiggles = [Wiggle(wig).parse() for wig in self.wiggle_files]
        base_columns_lst = []
        # Initialize matrix with full length genome
        for seq_record in parsed_fasta:
            for i in range(1, len(seq_record.seq) + 1, 1):
                base_columns_lst.append([seq_record.id, i])
        self.wiggle_matrix_df = pd.DataFrame(data=base_columns_lst, columns=["seqid", "location"])
        for parsed_wiggle in parsed_wiggles:
            condition_name = parsed_wiggle.at[0, "track_name"]
            self.wiggle_matrix_df = pd.merge(how='outer',
                                             left=self.wiggle_matrix_df,
                                             right=parsed_wiggle[["variableStep_chrom", "location", "score"]],
                                             left_on=['seqid', 'location'],
                                             right_on=['variableStep_chrom', 'location']).fillna(0.0)
            self.wiggle_matrix_df.rename(columns={"score": condition_name}, inplace=True)
            for column in self.wiggle_matrix_df.columns.tolist():
                if "variableStep_chrom" in column:
                    self.wiggle_matrix_df.drop(column, axis=1, inplace=True)
        return self.wiggle_matrix_df

    def get_matrix_by_orientation(self):
        self.build_matrix()
        f_column_list = ["seqid", "location"]
        r_column_list = ["seqid", "location"]
        for column in self.wiggle_matrix_df.columns.tolist():
            if "seqid" != column != "location":
                if self.wiggle_matrix_df[self.wiggle_matrix_df[column] < 0].empty:
                    f_column_list.append(column)
                if self.wiggle_matrix_df[self.wiggle_matrix_df[column] > 0].empty:
                    r_column_list.append(column)
        return self.wiggle_matrix_df[f_column_list], self.wiggle_matrix_df[r_column_list]