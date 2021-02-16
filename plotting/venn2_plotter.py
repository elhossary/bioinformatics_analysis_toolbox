from matplotlib_venn import venn2
from matplotlib import pyplot as plt

A = 126
B = 44
AB_intersect = 10
subsets = (A - AB_intersect, B - AB_intersect, AB_intersect)
labels = (f"{A} E. coli's binding sites upstream of genes + TSS\n(having protein clusters)",
          f"{B} K. pneumonia binding sites upstream of genes + TSS\n(having protein clusters)")
fig = plt.figure(figsize=(12, 4))
venn2(subsets=subsets, set_labels=labels, alpha=0.5)
plt.title("Genes belonging to the same protein clusters + TSS located downstream of binding sites")
fig.savefig("orth2.png")
