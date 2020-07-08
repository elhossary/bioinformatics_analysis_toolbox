from wiggle import Wiggle
import pandas as pd
import os
import glob
from Bio import SeqIO
import multiprocessing as mp


class WiggleMatrix:
    def __init__(self, fasta_file, wiggle_files, processes=None):
        self.fasta_file = fasta_file
        self.wiggle_files = wiggle_files
        self.wiggle_matrix_df = None
        if processes is None:
            self.processes = mp.cpu_count()
        elif processes.isnumeric():
            processes = int(processes)
            if processes > mp.cpu_count():
                self.processes = mp.cpu_count()
            elif processes > 0:
                self.processes = processes
            else:
                print("Error in processes number")
                exit(1)
        else:
            print("Error in processes number")
            exit(1)


    def build_matrix(self):
        parsed_fasta = SeqIO.parse(self.fasta_file, "fasta")
        pool = mp.Pool(processes=self.processes)
        processes = []
        for wig in self.wiggle_files:
            processes.append(pool.apply_async(self._parse_single_wiggle, (wig, )))
        parsed_wiggles = [p.get() for p in processes]
        pool.close()
        base_columns_lst = []
        # Initialize matrix with full length genome
        total_len = 0
        # Build dataframe backbone
        for seq_record in parsed_fasta:
            base_columns_lst = [[seq_record.id, i] for i in range(1, len(seq_record.seq) + 1, 1)]
            total_len += len(seq_record.seq)
        self.wiggle_matrix_df = pd.DataFrame(data=base_columns_lst, columns=["seqid", "location"])
        pool = mp.Pool(processes=self.processes)
        processes = []
        print("Building the matrix from the parsed files")
        for parsed_wiggle in parsed_wiggles:
            processes.append(pool.apply_async(self._merge_single_wiggle_to_matrix, (parsed_wiggle,)))
        columns_series = [p.get() for p in processes]
        for column in columns_series:
            self.wiggle_matrix_df[column.name] = column
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

    @staticmethod
    def _parse_single_wiggle(wig):
        return Wiggle(wig).parse()

    def _merge_single_wiggle_to_matrix(self, wig):
        condition_name = wig.at[0, "track_name"]
        self.wiggle_matrix_df = pd.merge(how='outer',
                                         left=self.wiggle_matrix_df,
                                         right=wig[["variableStep_chrom", "location", "score"]],
                                         left_on=['seqid', 'location'],
                                         right_on=['variableStep_chrom', 'location']).fillna(0.0)
        self.wiggle_matrix_df.rename(columns={"score": condition_name}, inplace=True)
        for column in self.wiggle_matrix_df.columns.tolist():
            if "variableStep_chrom" in column:
                self.wiggle_matrix_df.drop(column, axis=1, inplace=True)
        return self.wiggle_matrix_df[condition_name]

    def write_matrix_to_wiggle_files(self, matrix, out_dir, prefix=None):
        print("==> Writing wiggle files")
        seqids = matrix["seqid"].unique()
        columns = [col for col in matrix.columns if col not in ["seqid", "location"]]
        for col in columns:
            out_str = ""
            out_str += f'track type=wiggle_0 name="{col}"\n'
            for seqid in seqids:
                out_str += f'variableStep chrom={seqid} span=1\n'
                out_str += matrix[matrix["seqid"] == seqid][["location", col]]\
                    .to_csv(index=False, header=False, mode='a', sep=" ")
            file = open(os.path.abspath(f"{os.path.dirname(out_dir)}/{prefix}{col}.wig"), "w")
            file.write(out_str)
            file.close()
            s = '\n     └── '
            print(f"===> Wrote file: {prefix}{col}.wig, contains sequence IDs:\n"
              f"     └── {s.join(seqids)}")