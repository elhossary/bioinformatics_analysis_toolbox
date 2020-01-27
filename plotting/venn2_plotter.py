from matplotlib_venn import venn2
from matplotlib import pyplot as plt

subsets = (2501, 2449, 718)
# (2501, 2449, 718)
# A, B, AB intersect
labels = ("Terminators based on Term-Seq",
          "Terminators based on dRNA-Seq")
fig = plt.figure(figsize=(8, 4))
venn2(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("Terminators overlapping plot")
fig.savefig("term_venn.png")
