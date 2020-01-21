from wiggle_parser import WiggleParser as wp
import matplotlib.pyplot as plt
import pandas as pd
import collections
import numpy as np
from numpy import genfromtxt


def main():
	f_wigs_parsed = wp("test_data/cholerae/WT_S2_0.1_plus_TEX_forward.wig").parse()
	r_wigs_parsed = wp("test_data/cholerae/WT_S2_0.1_plus_TEX_forward.wig").parse()
	tss_arr = build_arr_form_gff("../test_data/cholerae/GCF_000006745.1_ASM674v1_genomic_TSS.gff")
	scores = []
	all = tss_arr.shape[0]
	not_found = 0
	for key, value in f_wigs_parsed.items():
		for tss_index, tss_row in enumerate(tss_arr):
			if key == tss_row[0]:
				if tss_row[6] == '+':
					try:
						scores.append(value.iat[value.index[value[0] == tss_row[3]].values[0], 1])
					except:
						not_found += 1
						pass
	for key, value in r_wigs_parsed.items():
		for tss_index, tss_row in enumerate(tss_arr):
			if key == tss_row[0]:
				if tss_row[6] == '-':
					try:
						scores.append(value.iat[value.index[value[0] == tss_row[3]].values[0], 1].abs())
					except:
						not_found += 1
						pass

	start_at = 0
	max_range = 200
	resolution = 5
	fig = plt.figure(figsize=(16, 9))
	plt.hist([i for i in scores if start_at < i <= max_range], bins=2000)
	plt.title(f"TSS positions coverage distribution from raw TEX+ - Total {all}, Found {len(scores)}, Missing {not_found}")
	plt.xlabel(f"Coverage range from {start_at} to {max_range}")
	plt.ylabel("Frequency")
	plt.xticks(range(start_at, max_range + 1, resolution))
	plt.grid(True)
	fig.savefig(f"{start_at}-{max_range}_TSS_positions_coverage_distribution.png")


def build_arr_form_gff(path):
	data_arr = genfromtxt(path, delimiter="\t", comments="#", dtype=None, encoding=None)
	return data_arr

main()