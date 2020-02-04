from matplotlib_venn import venn3
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--title", required=True, nargs='+', help="", type=str)
parser.add_argument("--A_all", required=True, help="", type=int)
parser.add_argument("--A_all_title", required=True, nargs='+', help="", type=str)
parser.add_argument("--B_all", required=True, help="", type=int)
parser.add_argument("--B_all_title", required=True, nargs='+', help="", type=str)
parser.add_argument("--C_all", required=True, help="", type=int)
parser.add_argument("--C_all_title", required=True, nargs='+', help="", type=str)
parser.add_argument("--AB", required=True, help="AB_intersect", type=int)
parser.add_argument("--AC", required=True, help="AC_intersect", type=int)
parser.add_argument("--BC", required=True, help="BC_intersect", type=int)
parser.add_argument("--ABC", required=True, help="ABC_intersect", type=int)
parser.add_argument("--output", required=True, help="Path to output file with extension PNG or PDF", type=str)
args = parser.parse_args()

args.title = ' '.join(args.title)
args.A_all_title = ' '.join(args.A_all_title)
args.B_all_title = ' '.join(args.B_all_title)
args.C_all_title = ' '.join(args.C_all_title)
ABnotC = args.AB - args.ABC
BCnotA = args.BC - args.ABC
ACnotB = args.AC - args.ABC
A = args.A_all - (args.ABC + ABnotC + ACnotB)
B = args.B_all - (args.ABC + ABnotC + BCnotA)
C = args.C_all - (args.ABC + BCnotA + ACnotB)
subsets = (A, B, ABnotC, C, ACnotB, BCnotA, args.ABC)
labels = (f"{args.A_all} {args.A_all_title}",
          f"{args.B_all} {args.B_all_title}",
          f"{args.C_all} {args.C_all_title}")
fig = plt.figure(figsize=(10, 5))
venn3(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title(f"{args.A_all + args.B_all + args.C_all} {args.title}")
fig.savefig(args.output)
