import pandas as pd
from io import StringIO


class GFF_Overlap_Merger:

    def __init__(self, gff_str, annotation_type, merge_range, annotate):
        self.gff_str = gff_str
        self.merge_range = merge_range
        self.annotation_type = annotation_type
        self.annotate = annotate

    def merge_overlaps(self):
        col_names = ["accession", "source", "type", "start", "end", "dot1", "strand", "dot2", "attributes"]
        ret_gff_str = ""
        gff_df = pd.read_csv(StringIO(self.gff_str), names=col_names, sep="\t", comment="#").dropna()
        accession_list = list(gff_df.accession.unique())
        df_dict = {}
        for acc in accession_list:
            df_dict[f"{acc}_f"] = \
                self.merge_interval_lists(gff_df[(gff_df['accession'] == acc) & (gff_df['strand'] == "+")]
                                          .loc[:, ['start', 'end']].sort_values(by=['start', 'end']).values.tolist(),
                                          self.merge_range, self.annotate)
            df_dict[f"{acc}_r"] = \
                self.merge_interval_lists(gff_df[(gff_df['accession'] == acc) & (gff_df['strand'] == "-")]
                                          .loc[:, ['start', 'end']].sort_values(by=['start', 'end']).values.tolist(),
                                          self.merge_range, self.annotate)
        seq_types = gff_df.type.unique().tolist()
        if self.annotation_type == "":
            if len(seq_types) == 1:
                self.annotation_type = seq_types[0]
            else:
                self.annotation_type = "unknown_sequence_type"
        strand_func = lambda x: "+" if "_f" in x else "-"
        strand_letter_func = lambda x: "F" if "+" in x else "R"
        for acc in accession_list:
            for dict_key in df_dict.keys():
                if dict_key == f"{acc}_f" or dict_key == f"{acc}_r":
                    for loc in df_dict[dict_key]:
                        ret_gff_str += \
                            f"{acc}\t" + \
                            f"GFF_merger\t" + \
                            f"merged_{self.annotation_type}_seq\t" + \
                            f"{loc[0]}\t" + \
                            f"{loc[1]}\t" + \
                            f".\t" + \
                            f"{strand_func(dict_key)}\t" + \
                            f".\t" + \
                            f".\n"
        ret_gff_df = pd.read_csv(StringIO(ret_gff_str), names=col_names, sep="\t", comment="#")
        ret_gff_df = ret_gff_df.sort_values(by=['accession', 'start'])
        output_counter = 0
        last_accession = ""
        # Writing attributes
        for index, row in ret_gff_df.iterrows():
            if last_accession != row['accession']:
                last_accession = row['accession']
                output_counter = 0
            output_counter += 1
            ret_gff_df.at[index, 'attributes'] = f"id={row['accession']}_{strand_letter_func(row['strand'])}"\
                                                 + f"_{self.annotation_type}_{output_counter};"\
                                                 + f"name={row['accession']}_{strand_letter_func(row['strand'])}"\
                                                 + f"_{self.annotation_type}_{output_counter};"\
                                                 + f"seq_len={row['end'] - row['start']}"

        ret_gff_str = ret_gff_df.to_csv(sep="\t", index=False, header=False)
        return ret_gff_str, gff_df.shape[0], ret_gff_df.shape[0]

    @staticmethod
    def merge_interval_lists(list_in, merge_range, annotate):
        merge_range += 2
        list_out = []
        overlap_indices = []
        for loc in list_in:
            if len(list_out) == 0:
                list_out.append(loc)
            else:
                if loc[0] in range(list_out[-1][0], list_out[-1][-1] + merge_range):
                    list_out[-1][-1] = max([loc[-1], list_out[-1][-1]])
                    overlap_indices.append(list_out.index(list_out[-1]))
                else:
                    list_out.append(loc)
        if annotate == 'overlaps':
            list_out = [list_out[i] for i in overlap_indices]
        if annotate == 'no_overlaps':
            for i in overlap_indices:
                del list_out[i]
        return list_out
