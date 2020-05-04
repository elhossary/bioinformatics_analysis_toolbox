import pandas as pd


class GffParser:

    def __init__(self, path):
        self.path = path

    def parse(self):
        col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        df_dict = {}
        gff_df = pd.read_csv(self.path, names=col_names, sep="\t", comment="#").dropna()
        seqid_list = list(gff_df.seqid.unique())
        for seqid in seqid_list:
            df_dict[seqid] = gff_df[(gff_df['seqid'] == seqid)]
        return df_dict
