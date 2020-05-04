from more_itertools import consecutive_groups
import pandas as pd
from numpy import diff, where, split
from Bio import SeqIO
import glob
import matplotlib.pyplot as plt
from os import path


class PolyTStretchFinder:
    def __init__(self, refseq_files):
        self.refseq_files = refseq_files

    def find_stretches(self, min_len, t_content):
        ret_df = pd.DataFrame(columns=['seqid', 'start', 'end', 'strand', 't_content'])
        for file in glob.glob(self.refseq_files):
            fasta_parsed = SeqIO.parse(file, "fasta")
            for seq_record in fasta_parsed:
                print(f"Finding poly-T stretches for sequence {seq_record.id}")
                ret_df = ret_df.append(self.group_positions(str(seq_record.seq), seq_record.id, min_len, t_content))
        return ret_df

    def group_positions(self, sequence_str, seqid, min_len, t_content):
        base = "T"
        complement_base = lambda x: "T" if base == "A" else "A"
        # Get all indexes for T
        f_indices = [i for i, a in enumerate(sequence_str, 1) if a == base]
        r_indices = [i for i, a in enumerate(sequence_str, 1) if a == complement_base(base)]
        f_indices.sort()
        r_indices.sort()
        # Group continuous occurrences of T indexes in list
        f_signals = [list(group) for group in consecutive_groups(f_indices)]
        r_signals = [list(group) for group in consecutive_groups(r_indices)]
        # Group interrupted occurrences of T indexes in list and append
        for max_interruption in range(1, 3, 1):
            f_signals.extend(list(map(list, split(f_indices, where(diff(f_indices) > max_interruption + 1)[0] + 1))))
            r_signals.extend(list(map(list, split(r_indices, where(diff(r_indices) > max_interruption + 1)[0] + 1))))


        # Drop invalid signals
        f_signals = [self.validate_stretch(sig, min_len, t_content) for sig in f_signals]
        #f_signals = [sig for sig in f_signals if sig is not None]
        r_signals = [self.validate_stretch(sig, min_len, t_content) for sig in r_signals]
        #r_signals = [sig for sig in r_signals if sig is not None]
        stretches = []
        for sig in f_signals:
            if sig is None:
                continue
            sig_len = sig[-1] - sig[0] + 1
            sig_t_content = round(len(sig) / sig_len, 2)
            stretch = [seqid, sig[0], sig[-1], "+", sig_t_content, int(sig_len)]
            if stretch not in stretches:
                stretches.append(stretch)
        for sig in r_signals:
            if sig is None:
                continue
            sig_len = sig[-1] - sig[0] + 1
            sig_t_content = round(len(sig) / sig_len, 2)
            stretch = [seqid, sig[0], sig[-1], "-", sig_t_content, int(sig_len)]
            if stretch not in stretches:
                stretches.append(stretch)

        df_header = ['seqid', 'start', 'end', 'strand', 't_content', 'length']
        return pd.DataFrame.from_records(stretches, columns=df_header).sort_values(['seqid', 'start', 'end'])

    def trim_start(self, lst, min_len):
        if lst is not None:
            if len(lst) < min_len:
                return None
            if lst[1] - lst[0] > 1:
                del lst[0]
                if len(lst) == lst[-1] - lst[0] == min_len:
                    return lst
                else:
                    return self.trim_start(lst, min_len)
        return lst

    def trim_end(self, lst, min_len):
        if lst is not None:
            if len(lst) < min_len:
                return None
            if lst[-1] - lst[-2] > 1:
                del lst[-1]
                if len(lst) == lst[-1] - lst[0] == min_len:
                    return lst
                else:
                    return self.trim_end(lst, min_len)
        return lst

    def validate_stretch(self, signal, min_len, t_content):
        if len(signal) == signal[-1] - signal[0] == min_len:
            pass
        else:
            signal = self.trim_start(signal, min_len)
            signal = self.trim_end(signal, min_len)
        # validate after trimming
        if signal is not None:
            sig_t_count = len(signal)
            sig_t_content = sig_t_count / (signal[-1] - signal[0] + 1)
            if sig_t_count < min_len or sig_t_content < t_content:
                return None
        return signal

    def write_to_gff(self, data_df):
        print("Writing GFF file...")
        str_out = ""
        count = 0
        strand_letter_func = lambda x: "F" if x == "+" else "R"
        for index, row in data_df.iterrows():
            count += 1
            str_out += \
                f"{row['seqid']}\t" + \
                f"Poly_T_stretch_finder\t" + \
                f"Poly_T_stretch\t" + \
                f"{row['start']}\t" + \
                f"{row['end']}\t" + \
                f".\t" + \
                f"{row['strand']}\t" + \
                f".\t" + \
                f"id={row['seqid']}_{strand_letter_func(row['strand'])}_poly_T_stretch_{count};" + \
                f"name={row['seqid']}_{strand_letter_func(row['strand'])}_poly_T_stretch_{count};" + \
                f"length={row['length']};t_content={row['t_content']};" + \
                "\n"
        outfile = open(f"{path.splitext(path.basename(self.refseq_files))[0]}_poly_t_stretches.gff", "w")
        outfile.write(f"###gff-version 3\n{str_out}###")
        outfile.close()

    def generate_stats(self, in_df, data_type):
        print(f"Generating {data_type} stats...")
        val_list = in_df[data_type].values.tolist()
        bins = len(list(set(val_list)))
        fig = plt.figure(figsize=(8, 6))
        plt.hist(val_list, bins=bins)
        plt.title(f"Poly-T stretches {data_type} frequencies, total: {len(val_list)}, average: {sum(val_list) / len(val_list)}")
        plt.xlabel(f"{data_type}")
        plt.ylabel("Frequency")
        plt.grid(True)
        fig.savefig(f"{path.splitext(path.basename(self.refseq_files))[0]}_poly_T_{data_type}_frequencies.png")
