from wiggle_parser.wiggle_parser import WiggleParser as wp
import matplotlib.pyplot as plt
import pandas as pd
import collections


def main():
	f_wigs_parsed = wp("test_data/cholerae/TERM-SEQ_LAST-BASE_12898-BCM-rep2_S38_R1_001_div_by_12099864.0_multi_by_10181443.0_forward.wig").parse()
	r_wigs_parsed = wp("test_data/cholerae/TERM-SEQ_LAST-BASE_12898-BCM-rep2_S38_R1_001_div_by_12099864.0_multi_by_10181443.0_reverse.wig").parse()
	all_arr = pd.DataFrame()
	for key, value in f_wigs_parsed.items():
		all_arr = all_arr.append(value, ignore_index=True)
	for key, value in r_wigs_parsed.items():
		all_arr = all_arr.append(value.abs(), ignore_index=True)
	start_at = 20
	max_range = 200
	resolution = 5
	filtered_list = all_arr[all_arr[1].between(start_at, max_range)][1].to_list()
	fig = plt.figure(figsize=(16, 9))
	plt.hist(filtered_list, bins=200)
	plt.title("Coverage distribution for normalized the Term-Seq last-base")
	plt.xlabel(f"Coverage range from {start_at} to {max_range}")
	plt.ylabel("Frequency")
	plt.xticks(range(start_at, max_range + 1, resolution))
	plt.grid(True)
	fig.savefig(f"{start_at}-{max_range}_coverage_distribution.png")

main()
