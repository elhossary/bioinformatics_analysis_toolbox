import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


def decode_set_name(set_name):
    set_list = set_name.split("_")
    r_str = ""
    for i in set_list:
        if "mi" in i:
            r_str += i.replace("mi", "Max interruption: ") + ", "
        if "ws" in i:
            r_str += i.replace("ws", "Window size: ") + ", "
        if "t" in i:
            r_str += i.replace("t", "Tolerance: ") + ", "
        if "ml" in i:
            r_str += i.replace("ml", "TSS-terminator distance: ")
    return r_str


df = pd.read_csv("V_cholerae_O1_N16961_param_sets_results.csv", sep="\t")
set_names = df.param_set_name.values.tolist()
set_names = list(set([re.sub('co.*?_', '', x, flags=re.DOTALL).replace("co_", "") for x in set_names]))
cutoffs = [float(x) for x in list(set(df.min_coverage.values.tolist()))]
labels = [f"Cutoff {x}" for x in cutoffs]

merged_srna = []
tp = []
fn = []
fp = []
golden_set = df.golden_set.values.tolist()[0]
for set_name in set_names:
    for cutoff in cutoffs:
        for index, row in df.iterrows():
            if f"co{cutoff}_{set_name}" == row['param_set_name']:
                merged_srna.append(row['merged_srna_count'])
                tp.append(row['TP'])
                fn.append(row['FN'])
                fp.append(row['FP'])

    x = np.arange(len(labels))  # the label locations
    width = 0.2  # the width of the bars

    fig, ax = plt.subplots(figsize=(8, 6))
    rects_merged_srna = ax.bar(x + width, merged_srna, width, label='Merged sRNAs')
    rects_tp = ax.bar(x + width * 2, tp, width, label='True Positives')
    rects_fn = ax.bar(x + width * 3, fn, width, label='False Negatives')
    rects_fp = ax.bar(x + width * 4, fp, width, label='False Positives')
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Counts')
    title=f"sRNA predictions comparison to golden set ({golden_set})\n" \
          f"Variable parameter: minimum coverage cutoff\n" \
          f"Fixed parameters variable:\n{decode_set_name(set_name)}"

    ax.set_title(title)
    ax.set_xticks(x + width * 2.50)
    ax.set_xticklabels(labels)
    ax.legend()
    plt.yticks(list(range(0, max(fp) + 10, 50)))
    autolabel(rects_merged_srna)
    autolabel(rects_tp)
    autolabel(rects_fn)
    autolabel(rects_fp)
    fig.tight_layout()
    plt.savefig(f"count_stats_plots/count_stats_{set_name}_compare.png")
    plt.close(fig)
    merged_srna = []
    tp = []
    fn = []
    fp = []