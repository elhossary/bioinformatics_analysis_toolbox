from matplotlib_venn import venn2
from matplotlib import pyplot as plt

A = 3089
B = 1340
AB_intersect = 264
subsets = (A - AB_intersect, B - AB_intersect, AB_intersect)
labels = (f"{A} Terminators\nbased on dRNA-Seq", f"{B} Terminators\nbased on Term-Seq")
fig = plt.figure(figsize=(8, 4))
venn2(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("Terminators overlapping plot")
fig.savefig("ecoli_term_venn_cutoff_11.png")
