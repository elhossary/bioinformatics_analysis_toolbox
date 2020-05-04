from more_itertools import consecutive_groups
import pandas as pd

class PolyBaseStretchFinder:
    def __init__(self, forward_sequence_str, base_on_forward, orientation, tolerance=0):
        self.forward_sequence_str = forward_sequence_str
        self.base_on_forward = base_on_forward
        self.orientation = orientation
        self.tolerance = tolerance
        self.signals_out = []

    def get_signals(self):
        if self.tolerance > 0:
            self.signals_out = self.merge_interval_lists(self.get_consecutive_bases())
        elif self.tolerance == 0:
            self.signals_out = self.get_consecutive_bases()
        else:
            print("Error: merge range can't be a negative value")

    def get_consecutive_bases(self, seq_str, center_loc, range, base, min_len=5):
        #seq_str = seq_str[center_loc - range: center_loc + range]

        indices = [i for i, a in enumerate(seq_str, center_loc - range + 1) if a == base]
        indices.sort()
        all_signals = [list(group) for group in consecutive_groups(indices)]
        all_signals = [[i[0], i[-1]] for i in all_signals if len(i) >= min_len]
        return all_signals

    def merge_interval_lists(self, list_in, merge_range):
        merge_range += 2
        list_out = []
        for loc in list_in:
            if len(list_out) == 0:
                list_out.append(loc)
            else:
                if loc[0] in range(list_out[-1][0], list_out[-1][-1] + merge_range):
                    list_out[-1][-1] = loc[-1]
                else:
                    list_out.append(loc)
        return list_out

    def write_to_gff(self, seqid, output_file):
        print("Writing GFF file...")
        term_gff_str = ""
        count = 0
        current_accession = ""
        out_df = pd.DataFrame.from_records(self.signals_out)
        out_df = out_df.sort_values([0, 1])
        strand_letter_func = lambda x: "F" if x == "+" else "R"
        for index, row in out_df.iterrows():
            if current_accession != row[0] or current_accession == "":
                # term_gff_str += "###\n"
                current_accession = row[0]
                count = 0
            count += 1
            term_gff_str += \
                f"{row[0]}\t" + \
                f"Poly_T_stretch_finder\t" + \
                f"Poly_T_stretch\t" + \
                f"{row[1]}\t" + \
                f"{row[2]}\t" + \
                f".\t" + \
                f"{row[3]}\t" + \
                f".\t" + \
                f"id={current_accession}_{strand_letter_func(row[3])}_poly_T_terminator_{count};" + \
                f"name={current_accession}_{strand_letter_func(row[3])}_poly_T_terminator_{count}\n"
        outfile = open(args.gff_out, "w")
        outfile.write(f"###gff-version 3\n{term_gff_str}###")
        outfile.close()