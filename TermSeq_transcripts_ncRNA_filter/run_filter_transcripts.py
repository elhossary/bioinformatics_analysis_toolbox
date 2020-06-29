from Bio import SeqIO
from numpy import diff, where, split
import argparse
import pandas as pd
import glob
from gff_parser import GffParser
from poly_t_stretch_finder import PolyTStretchFinder as t_finder


def main():
    # Params
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_in", required=True, help="RefSeq fasta file", type=str)
    parser.add_argument("--termseq_transcripts_gff", required=True, help="", type=str)
    parser.add_argument("--range", default=20, required=False, help="", type=int)
    parser.add_argument("--min_len", default=20, required=False, help="", type=int)
    parser.add_argument("--max_len", default=300, required=False, help="", type=int)

    args = parser.parse_args()
    # ---------------------------
    """
    t_finder_obj = t_finder(args.fasta_in)
    stretches_df = t_finder_obj.find_stretches(5, 0.7)
    t_finder_obj.write_to_gff(stretches_df)
    t_finder_obj.generate_stats(stretches_df, "length")
    t_finder_obj.generate_stats(stretches_df, "t_content")
    """
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    #fasta_parsed = SeqIO.parse(glob.glob(args.fasta_in)[0], "fasta")

    gff_parsed = pd.read_csv(args.termseq_transcripts_gff, names=col_names, sep="\t", comment="#")
    stretches_df = pd.read_csv("GCF_000006745.1_ASM674v1_genomic_term_TransTermHP_only.gff", names=col_names, sep="\t", comment="#")

    output_df = pd.DataFrame()
    counter = 0
    for index, transcript in gff_parsed.iterrows():
        if args.min_len <= int(transcript['end']) - int(transcript['start']) + 1 <= args.max_len:
            counter += 1
        if transcript['strand'] == "+":
            selection = stretches_df[(stretches_df['seqid'] == transcript['seqid']) &
                                     (stretches_df['strand'] == transcript['strand']) &
                                     (stretches_df['end'] >= transcript['end'] - args.range) &
                                     (stretches_df['start'] <= transcript['end'] + args.range)]
            if not selection.empty:
                for i, selection_item in selection.iterrows():
                    transcript['end'] = selection_item['end']
                    output_df = output_df.append(transcript, ignore_index=True)
        elif transcript['strand'] == "-":
            selection = stretches_df[(stretches_df['seqid'] == transcript['seqid']) &
                                     (stretches_df['strand'] == transcript['strand']) &
                                     (stretches_df['end'] >= transcript['start'] - args.range) &
                                     (stretches_df['start'] <= transcript['start'] + args.range)]

            if not selection.empty:
                for i, selection_item in selection.iterrows():
                    transcript['start'] = selection_item['start']
                    output_df = output_df.append(transcript, ignore_index=True)
        else:
            print("Fatal error")
    print(counter)
    output_df.reset_index(inplace=True)
    write_to_gff(output_df, "Length_TransTermHP_filtered_TermSeq_transcript.gff", args.min_len, args.max_len)


def write_to_gff(data_df, gff_out, min_len=20, max_len=300):
    print("Writing GFF file...")
    str_out = ""
    count = 0
    strand_letter_func = lambda x: "F" if x == "+" else "R"
    for index, row in data_df.iterrows():
        if min_len <= int(row['end']) - int(row['start']) + 1 <= max_len:
            count += 1
            str_out += \
                f"{row['seqid']}\t" + \
                f"Filtered_TermSeq_Transcript\t" + \
                f"TermSeq_transcript\t" + \
                f"{int(row['start'])}\t" + \
                f"{int(row['end'])}\t" + \
                f".\t" + \
                f"{row['strand']}\t" + \
                f".\t" + \
                f"id={row['seqid']}_{strand_letter_func(row['strand'])}_termseq_transcript_{count};" + \
                f"name={row['seqid']}_{strand_letter_func(row['strand'])}_termseq_transcript_{count};" + \
                "\n"
    outfile = open(gff_out, "w")
    outfile.write(f"###gff-version 3\n{str_out}###")
    outfile.close()
    print(count)


main()