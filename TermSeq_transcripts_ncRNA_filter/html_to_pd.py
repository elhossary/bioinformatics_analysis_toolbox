import pandas as pd
import numpy as np

def getindices(s):
    return [i for i, c in enumerate(s) if c.isupper()]

fpath = "/home/muhoss/RegulonDB.html"

df_list = pd.read_html(fpath)

df = pd.DataFrame()
for i in df_list:
    df = pd.concat([df, i])

df = df.dropna(axis=0, subset=["Promoter Sequence"])
outstr1 = ""
outstr2 = ""
for idx in df.index:
    if "[1]" in df.at[idx, 'Reference(s)']:
        seq = str(df.at[idx, 'Promoter Sequence']).split(' ')[0]
        up_idx = getindices(seq)[0]
        if seq[:up_idx].upper() not in outstr1:
            outstr1 += f">{df.at[idx, 'Promoter Name']}\n" \
                      f"{seq[:up_idx].upper()}\n"
    else:
        seq = str(df.at[idx, 'Promoter Sequence']).split(' ')[0]
        up_idx = getindices(seq)[0]
        if seq[:up_idx].upper() not in outstr2:
            outstr2 += f">{df.at[idx, 'Promoter Name']}\n" \
                      f"{seq[:up_idx].upper()}\n"
with open("/home/muhoss/SigmaE_promoter_sequences_RegulonDB_Rhodius_only.fa", "w") as f:
    f.write(outstr1[:-1])
with open("/home/muhoss/SigmaE_promoter_sequences_RegulonDB_without_Rhodius.fa", "w") as f:
    f.write(outstr2[:-1])
with open("/home/muhoss/SigmaE_promoter_sequences_RegulonDB_all.fa", "w") as f:
    f.write(f"{outstr1}{outstr2[:-1]}")