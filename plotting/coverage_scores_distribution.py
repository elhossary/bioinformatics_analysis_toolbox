from wiggle_parser import WiggleParser as wp
import matplotlib.pyplot as plt
import pandas as pd
import collections
import numpy as np

def main():
	f_wigs_parsed = wp("test_data/cholerae/TERM-SEQ_LAST-BASE_12898-BCM-rep2_S38_R1_001_div_by_12099864.0_multi_by_10181443.0_forward.wig").parse()
	r_wigs_parsed = wp("test_data/cholerae/TERM-SEQ_LAST-BASE_12898-BCM-rep2_S38_R1_001_div_by_12099864.0_multi_by_10181443.0_reverse.wig").parse()
	all_arr = []
	ranges = []
	counts = []

	for key, value in f_wigs_parsed.items():
		all_arr.append(value[1].to_numpy())
	for key, value in r_wigs_parsed.items():
		value[1] = value[1].abs()
		all_arr.append(value[1].to_numpy())
	n_arr = np.concatenate(all_arr)
	#n_arr = n_arr[n_arr > 9]
	last = 0
	max_score = int(np.amax(n_arr)) + 1
	max_range = 50
	interval = 5
	for i in range(0, max_range, interval):
		counts.append(((last < n_arr) & (n_arr <= i)).sum())
		ranges.append(i)
		last = i

	#ranges.append(last + 5)
	#counts.append(((last < n_arr) & (n_arr <= max_score)).sum())
	max_count = max(counts)
	fig = plt.figure(figsize=(16, 9))
	plt.plot(ranges, counts)
	plt.plot(ranges, counts, ".")
	plt.title("Score distribution for the Term-Seq last-base normalized coverage")
	plt.xlabel(f"Score ranges, Ranges from 0 to {max_range}, interval {interval}")
	plt.ylabel("Frequency of scores in ranges")
	plt.xticks(range(0, max_range + 1, interval))
	plt.yticks(range(0, max_count, round(int(max_count/50), -3)))
	plt.grid(True)
	fig.savefig(f"0-{max_range}_coverage_score_distribution.png")

main()