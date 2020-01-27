from matplotlib_venn import venn3
from matplotlib import pyplot as plt

subsets = (503, 317, 54, 105, 36, 27, 107)
# (503, 317, 54, 105, 36, 27, 107)
# A, B, AB intersect, C, AC intersect, BC intersect, ABC intersect
labels = ("Term-Seq sRNAs", "dRNA-Seq sRNAs", "Kai's annotations")
venn3(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("Overlapping plot")
plt.savefig("venn.png")