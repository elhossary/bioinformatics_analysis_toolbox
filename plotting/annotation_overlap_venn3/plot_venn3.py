from matplotlib_venn import venn3
from matplotlib import pyplot as plt
import argparse
import pandas as pd
import os
from gff_overlap_merger import GFF_Overlap_Merger as gff_mrg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--title", required=True, nargs='+', help="", type=str)
    parser.add_argument("--gff_a", required=True, help="", type=str)
    parser.add_argument("--gff_a_label", required=True, nargs='+', help="", type=str)
    parser.add_argument("--gff_b", required=True, help="", type=str)
    parser.add_argument("--gff_b_label", required=True, nargs='+', help="", type=str)
    parser.add_argument("--gff_c", required=True, help="", type=str)
    parser.add_argument("--gff_c_label", required=True, nargs='+', help="", type=str)
    parser.add_argument("--output", required=True, help="Path to output file with extension PNG or PDF", type=str)
    args = parser.parse_args()


    # Pre-merges for getting unique sets
    pre_merged_str_a, _, A_all = \
        gff_mrg(open(os.path.abspath(args.gff_a), "r").read(), "premerged", 0, "all").merge_overlaps()
    pre_merged_str_b, _, B_all = \
        gff_mrg(open(os.path.abspath(args.gff_b), "r").read(), "premerged", 0, "all").merge_overlaps()
    pre_merged_str_c, _, C_all = \
        gff_mrg(open(os.path.abspath(args.gff_c), "r").read(), "premerged", 0, "all").merge_overlaps()

    # Overlap merges
    AB_merge_str = pre_merged_str_a + pre_merged_str_b
    AB_merge_str, _, AB = gff_mrg(AB_merge_str, "ab_merge", 0, "overlaps").merge_overlaps()
    BC_merge_str = pre_merged_str_b + pre_merged_str_c
    BC_merge_str, _, BC = gff_mrg(BC_merge_str, "bc_merge", 0, "overlaps").merge_overlaps()
    AC_merge_str = pre_merged_str_a + pre_merged_str_c
    AC_merge_str, _, AC = gff_mrg(AC_merge_str, "ac_merge", 0, "overlaps").merge_overlaps()
    ABC_merge_str = AB_merge_str + BC_merge_str + AC_merge_str
    ABC_merge_str, _, ABC = gff_mrg(ABC_merge_str, "abc_merge", 0, "overlaps").merge_overlaps()

    # Plot
    plot(args.title, A_all, B_all, C_all, AB, BC, AC, ABC, args.gff_a_label.replace('_', ' '),
         args.gff_b_label.replace('_', ' '), args.gff_c_label.replace('_', ' '), args.output)


def plot(title, A_all, B_all, C_all, AB, BC, AC, ABC, A_all_title, B_all_title, C_all_title, output):
    title = ' '.join(title)
    ABnotC = AB - ABC
    BCnotA = BC - ABC
    ACnotB = AC - ABC

    A = A_all - (ABC + ABnotC + ACnotB)
    B = B_all - (ABC + ABnotC + BCnotA)
    C = C_all - (ABC + BCnotA + ACnotB)
    subsets = (A, B, ABnotC, C, ACnotB, BCnotA, ABC)
    labels = (f"{A_all} {labels_wrapper(' '.join(A_all_title))}",
              f"{B_all} {labels_wrapper(' '.join(B_all_title))}",
              f"{C_all} {' '.join(C_all_title)}")
    fig = plt.figure(figsize=(10, 6))
    venn3(subsets=subsets, set_labels=labels, alpha=0.5)
    plt.title(f"{A_all + B_all + C_all} {title.replace('_', ' ')}")
    fig.savefig(output)


def labels_wrapper(label_str):
    spaces_indices = [i for i, a in enumerate(label_str) if a == " "]
    for i in range(4, len(spaces_indices), 5):
        label_str = label_str[:spaces_indices[i]] + "\n" + label_str[spaces_indices[i]+1:]
    return label_str

main()
