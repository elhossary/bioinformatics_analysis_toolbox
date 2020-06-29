import pandas as pd

df = pd.read_csv("combine_tracks.csv", sep="\t")
forward_colors = ["green", "red"]
reverse_colors = ["blue", "magenta"]
bw_col_names = df.columns
bw_col_names = [x for x in bw_col_names if "bigwig_name_" in x]
if len(bw_col_names) == len(forward_colors) == len(reverse_colors):
    out_str = ""
    for index, row in df.iterrows():
        out_str += f"[tracks.{row['track_name']}]\n"
        out_str += f"key={row['track_name']}\n"
        out_str += "type=MultiBigWig/View/Track/MultiWiggle/MultiXYPlot\n"
        out_str += "storeClass=MultiBigWig/Store/SeqFeature/MultiBigWig\n"
        out_str += "autoscale=local\n"
        out_str += "nonCont=true\n"
        if row['strand'] == 'f':
            out_str += f"description={row['track_name']}, {forward_colors[0]} for TEX-, and {forward_colors[1]} for TEX+\n"
            for bw_name in bw_col_names:
                out_str += f'urlTemplates += json:{{"url": "raw_data/bigwig/{row[bw_name]}.bw", "color": "{forward_colors[bw_col_names.index(bw_name)]}", "category" : "Combined TEX +/- tracks"}}\n'
        if row['strand'] == 'r':
            out_str += f"description={row['track_name']}, {reverse_colors[0]} for TEX-, and {reverse_colors[1]} for TEX+\n"
            for bw_name in bw_col_names:
                out_str += f'urlTemplates += json:{{"url": "raw_data/bigwig/{row[bw_name]}.bw", "color": "{reverse_colors[bw_col_names.index(bw_name)]}", "category" : "Combined TEX +/- tracks"}}\n'
file_out = open("tracks.conf", "w")
file_out.write(out_str)
file_out.close()
