import pandas as pd
from io import StringIO


class GFF_Overlap_Merger:

    def __init__(self, gff_str, annotation_type, merge_range, annotate):
        self.gff_str = gff_str
        self.merge_range = merge_range
        self.annotation_type = annotation_type
        self.annotate = annotate

    def merge_overlaps(self):
        col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        gff_df = pd.read_csv(StringIO(self.gff_str), names=col_names, sep="\t", comment="#").dropna()
        seqid_list = list(gff_df.seqid.unique())
        df_dict_list = []
        f_no_overlap = []
        r_no_overlap = []
        for seqid in seqid_list:
            seqid_gff_df_slice = gff_df[gff_df['seqid'] == seqid]
            f_seqid_gff_df_slice = seqid_gff_df_slice[seqid_gff_df_slice['strand'] == "+"]
            r_seqid_gff_df_slice = seqid_gff_df_slice[seqid_gff_df_slice['strand'] == "-"]
            f_list_in = f_seqid_gff_df_slice.loc[:, ['start', 'end']].sort_values(by=['start', 'end']).values.tolist()
            r_list_in = r_seqid_gff_df_slice.loc[:, ['start', 'end']].sort_values(by=['start', 'end']).values.tolist()
            counter = 0
            for iteration in range(0, self.merge_range + 1, 1):
                f_tmp_overlaps, f_no_overlap = self.merge_interval_lists(f_list_in, iteration)
                r_tmp_overlaps, r_no_overlap = self.merge_interval_lists(r_list_in, iteration)
                f_list_in = f_no_overlap
                r_list_in = r_no_overlap
                for x in f_tmp_overlaps:
                    counter += 1
                    df_dict_list.append({"seqid": seqid, "start": x[0], "end": x[1], "strand": "+",
                                         "merge_range": iteration, "id": counter})
                for x in r_tmp_overlaps:
                    counter += 1
                    df_dict_list.append({"seqid": seqid, "start": x[0], "end": x[1], "strand": "-",
                                         "merge_range": iteration, "id": counter})
            for x in f_no_overlap:
                counter += 1
                df_dict_list.append({"seqid": seqid, "start": x[0], "end": x[1], "strand": "+",
                                     "merge_range": "unmerged", "id": counter})
            for x in r_no_overlap:
                counter += 1
                df_dict_list.append({"seqid": seqid, "start": x[0], "end": x[1], "strand": "-",
                                     "merge_range": "unmerged", "id": counter})
        seq_types = gff_df.type.unique().tolist()
        if self.annotation_type == "":
            if len(seq_types) == 1:
                self.annotation_type = seq_types[0]
            else:
                self.annotation_type = "unknown_sequence_type"
        ret_gff_df = pd.DataFrame.from_records(df_dict_list, columns=col_names + ["merge_range", "id"])
        ret_gff_df["source"], ret_gff_df["type"] = "GFF_merger", f"merged_{self.annotation_type}_seq"
        ret_gff_df["score"], ret_gff_df["phase"] = ".", "."
        merge_stats_dict = ret_gff_df["merge_range"].value_counts().to_dict()
        merge_stats = ""
        total = 0
        for k in merge_stats_dict.keys():
            if k == "unmerged":
                merge_stats += f"Total unmerged annotations: {merge_stats_dict[k]}\n"
            else:
                merge_stats += f"Merged annotations at {k} range: {merge_stats_dict[k]}\n"
                total += merge_stats_dict[k]
        merge_stats += f"Total merged annotations: {total}"
        if self.annotate == "all":
            pass
        elif self.annotate == "overlaps":
            ret_gff_df = ret_gff_df[ret_gff_df["merge_range"] != "unmerged"]
        elif self.annotate == "no_overlaps":
            ret_gff_df = ret_gff_df[ret_gff_df["merge_range"] == "unmerged"]
        else:
            print("Fatal error")
            exit(1)
        ret_gff_df = ret_gff_df.apply(lambda row: self.generate_attributes(row), axis=1)
        ret_gff_df.drop(["merge_range", "id"], inplace=True, axis=1)
        ret_gff_df.sort_values(by=['seqid', 'start'])
        ret_gff_str = ret_gff_df.to_csv(sep="\t", index=False, header=False)
        return ret_gff_str, gff_df.shape[0], ret_gff_df.shape[0], merge_stats

    @staticmethod
    def merge_interval_lists(list_in, merge_range):
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
        overlap_indices = list(set(overlap_indices))
        overlap_indices.sort()
        overlaps_list_out = [list_out[i] for i in overlap_indices]
        return overlaps_list_out, [i for i in list_out if i not in overlaps_list_out]

    def generate_attributes(self, row):
        strand_letter_func = lambda s: "F" if "+" in s else "R"
        row["attributes"] = f"id={row['seqid']}_{strand_letter_func(row['strand'])}" \
                            f"_{self.annotation_type}_{row['id']}" \
                            f";name={row['seqid']}_{strand_letter_func(row['strand'])}" \
                            f"_{self.annotation_type}_{row['id']}" \
                            f";seq_len={row['end'] - row['start'] + 1}" \
                            f";merge_range={row['merge_range']}"
        return row