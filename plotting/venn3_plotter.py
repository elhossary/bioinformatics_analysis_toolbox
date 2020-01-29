from matplotlib_venn import venn3
from matplotlib import pyplot as plt


A_all = 317 # dRNA srna
B_all = 854 # Term-Seq srna
C_all = 105 # Kai srna
AB_intersect = 77
AC_intersect = 36
BC_intersect = 34
random_all_intersects = 133
sum_all_2_node_intersect = AB_intersect + AC_intersect + BC_intersect
ABC_intersect = sum_all_2_node_intersect - random_all_intersects

ABnotC = AB_intersect - ABC_intersect
BCnotA = BC_intersect - ABC_intersect
ACnotB = AC_intersect - ABC_intersect

A = A_all - AC_intersect - ABnotC
B = B_all - BC_intersect - ABnotC
C = C_all - BC_intersect - ACnotB

subsets = (A, B, ABnotC, C, ACnotB, BCnotA, ABC_intersect)


labels = (f"{A_all} sRNAs based on dRNA-Seq",
          f"{B_all} sRNAs based on dRNA-seq TSS\nand Term-Seq Terminators",
          f"{C_all} sRNAs from Kai's annotations")
fig = plt.figure(figsize=(10, 5))
venn3(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("sRNAs overlapping plot")
fig.savefig("sRNA_venn.png")
