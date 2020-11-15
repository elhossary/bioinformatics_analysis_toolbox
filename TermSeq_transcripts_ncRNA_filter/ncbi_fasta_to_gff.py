fpath = "/home/muhoss/Downloads/all_acc.fa"

gff_str = ""
gff_dict_list = []
lst = []
with open(fpath, "r") as f_content_file:
    for line in f_content_file.readlines():
        if line.startswith(">"):
            lst.append(line.split("|ref|")[1].split("|"))
lst = [[i[0], i[1].replace(":", "").replace("-", " ").split(" ")] for i in lst]
lst = [[i[0]] + i[1] for i in lst]
lst = [i[0:4] for i in lst]
for i in lst:
    i[1] = int(i[1].replace("c", ""))
    i[2] = int(i[2])
    if i[1] > i[2]:
        start = i[2]
        end = i[1]
        strand = "-"
    else:
        start = i[1]
        end = i[2]
        strand = "+"
    gff_dict_list.append(
        {
            "seqid" : i[0],
            "source": "NCBI",
            "type": "ORF",
            "start": start,
            "end": end,
            "score": ".",
            "strand": strand,
            "phase": ".",
            "attributes": f"ID={i[3]}{strand};Name={i[3]}{strand}"
        }
    )
import pandas
df = pandas.DataFrame(gff_dict_list)
df.sort_values(["seqid", "start", "end"], inplace=True, axis=0)

for seqid in df["seqid"].unique():
    df[df["seqid"] == seqid].to_csv(f"/home/muhoss/Downloads/{seqid}_ORFs.gff", sep="\t", header=False, index=False)
