import pandas as pd
df = pd.read_csv("V_cholerae_O1_N16961_param_sets_results.csv", sep="\t")
x = df.to_markdown()
y = open("table.md", "w")
y.write(x)
y.close()