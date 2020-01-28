from matplotlib_venn import venn2
from matplotlib import pyplot as plt

A = 2492
B = 2441
AB_intersect = 722
subsets = (A - AB_intersect, B - AB_intersect, AB_intersect)
labels = (f"{A} Terminators based on dRNA-Seq", f"{B} Terminators based on Term-Seq")
fig = plt.figure(figsize=(8, 4))
venn2(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("Terminators overlapping plot")
fig.savefig("term_venn.png")
