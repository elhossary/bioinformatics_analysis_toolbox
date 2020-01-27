import matplotlib.pyplot as plt
import numpy as np


labels = ['Paper published', 'TSSpredator default', 'ANNOgesic optimized']
tss = [7202, 12692, 2286]
transcripts = [8230, 9576, 8230]
terminators = [2501, 2471, 2501]
srnas = [317, 422, 464]

x = np.arange(len(labels))  # the label locations
width = 0.2  # the width of the bars

fig, ax = plt.subplots(figsize=(8, 6))
rects_tss = ax.bar(x + width, tss, width, label='TSS')
rects_transcripts = ax.bar(x + width * 2, transcripts, width, label='Transcripts')
rects_terminators = ax.bar(x + width * 3, terminators, width, label='Terminator')
rects_srnas = ax.bar(x + width * 4, srnas, width, label='sRNAs')
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Counts')
ax.set_title('ANNOgesic prediction counts for each parameter set')
ax.set_xticks(x + width * 2.50)
ax.set_xticklabels(labels)
ax.legend()
plt.yticks(list(range(0, max(tss) + 100, 500)))


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


autolabel(rects_tss)
autolabel(rects_transcripts)
autolabel(rects_terminators)
autolabel(rects_srnas)
fig.tight_layout()

plt.savefig("TEX_bar_chart.png")
