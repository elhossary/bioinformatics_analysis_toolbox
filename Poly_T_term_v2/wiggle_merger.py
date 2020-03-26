from wiggle_parser import WiggleParser as wp
import pandas as pd


class WiggleMerger:

    def __init__(self, wig_files, fasta_parsed):
        self.wig_files = wig_files
        self.fasta_parsed = fasta_parsed

    def merge_wigs_by_max(self):
        # Handling multiple wig files
        f_wigs_parsed = {}
        r_wigs_parsed = {}
        for seq_record in self.fasta_parsed:
            f_wigs_parsed[seq_record.id] = pd.DataFrame(data=range(1, len(seq_record.seq), 1))
            r_wigs_parsed[seq_record.id] = pd.DataFrame(data=range(1, len(seq_record.seq), 1))

        for wig in self.wig_files:
            parsed_wig = wp(wig).parse()
            for accession, coverage in parsed_wig.items():
                if accession in r_wigs_parsed.keys():
                    if coverage[coverage[1] < 0].empty:
                        f_wigs_parsed[accession] = pd.merge(how='outer', left=f_wigs_parsed[accession],
                                                            right=coverage,
                                                            left_on=0, right_on=0).fillna(0.0)
                    if coverage[coverage[1] > 0].empty:
                        r_wigs_parsed[accession] = pd.merge(how='outer', left=r_wigs_parsed[accession],
                                                            right=coverage,
                                                            left_on=0, right_on=0).fillna(0.0)
        for accession in f_wigs_parsed.keys():
            f_wigs_parsed[accession][1] = f_wigs_parsed[accession].iloc[:, 1:-1].max(axis=1)
            f_wigs_parsed[accession] = f_wigs_parsed[accession].iloc[:, [0, -1]]
            #f_wigs_parsed[accession] = f_wigs_parsed[accession][f_wigs_parsed[accession][1] != 0.0]
        for accession in r_wigs_parsed.keys():
            r_wigs_parsed[accession][1] = r_wigs_parsed[accession].iloc[:, 1:-1].min(axis=1)
            r_wigs_parsed[accession] = r_wigs_parsed[accession].iloc[:, [0, -1]]
            #r_wigs_parsed[accession] = r_wigs_parsed[accession][r_wigs_parsed[accession][1] != 0.0]
        return f_wigs_parsed, r_wigs_parsed