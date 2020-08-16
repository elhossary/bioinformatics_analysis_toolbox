import matplotlib.pyplot as plt
#data = {90:910, 91:777, 92:667, 93:565, 94:459, 95:364, 96:285, 97:196, 98:120 ,99:78}
data = {90:2871, 91:2395, 92:2054, 93:1830, 94:1594, 95:1354, 96:1076, 97:844, 98:593 ,99:350}
fig = plt.figure(figsize=(16, 9))
plt.bar(data.keys(), data.values())
plt.title("Counts of predicted transcripts per nth percentile normalized coverage in E. coli")
plt.xlabel("nth percentile")
plt.ylabel("Transcripts count")
plt.xticks(range(90, 100, 1))
plt.yticks(range(0, 3000, 100))
plt.savefig(f"percentile_predictions.png")
plt.close(fig)