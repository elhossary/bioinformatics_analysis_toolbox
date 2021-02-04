from matplotlib_venn import venn2
from matplotlib import pyplot as plt

A = 246
B = 148
AB_intersect = 62
subsets = (A - AB_intersect, B - AB_intersect, AB_intersect)
labels = (f"{A} E. coli's genes\nwith orthologs", f"{B} K. pneumonia genes\nwith orthologs")
fig = plt.figure(figsize=(8, 4))
venn2(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("Orthologs with SigE binding sites shared between E. coli in K. pneumonia")
fig.savefig("orth2.png")
