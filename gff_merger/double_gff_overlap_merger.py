import argparse
import os
import glob
import pandas as pd
from numpy import genfromtxt
from io import StringIO

class Double_GFF_Overlap_Merger:

    def __init__(self, gff1, gff2, annotation_type):
        self.gff1 = gff1
        self.gff2 = gff2
        self.annotation_type = annotation_type

    def merge_overlaps(self):
        col_names = ["accession", "source", "type", "start", "end", "dot1", "strand", "dot2", "attributes"]
        gff1_df = pd.read_csv(self.gff1, names=col_names, sep="\t", comment="#")
        gff2_df = pd.read_csv(self.gff2, names=col_names, sep="\t", comment="#")
        

    def build_arr_form_gff(path):
        data_arr = genfromtxt(path, delimiter="\t", comments="#", dtype=None, encoding=None)
        return data_arr