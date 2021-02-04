from matplotlib_venn import venn2
from matplotlib import pyplot as plt

A = 2686
B = 5302
AB_intersect = 2083
subsets = (A - AB_intersect, B - AB_intersect, AB_intersect)
labels = (f"{A} E. coli's CDS with cluster IDs", f"{B} K. pneumonia CDS with cluster IDs")
fig = plt.figure(figsize=(8, 4))
venn2(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("Orthologs of E. coli in K. pneumonia")
fig.savefig("orth.png")
