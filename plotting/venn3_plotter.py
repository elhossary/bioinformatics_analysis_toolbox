from matplotlib_venn import venn3
from matplotlib import pyplot as plt

subsets = (503, 317, 54, 105, 36, 27, 107)
# (503, 317, 54, 105, 36, 27, 107)
# A, B, AB intersect, C, AC intersect, BC intersect, ABC intersect
labels = ("sRNAs based on dRNA-seq TSS\nand Term-Seq Terminators",
          "sRNAs based on dRNA-Seq",
          "sRNAs from Kai's annotations")
fig = plt.figure(figsize=(9, 6))
venn3(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("sRNAs overlapping plot")
fig.savefig("sRNA_venn.png")
