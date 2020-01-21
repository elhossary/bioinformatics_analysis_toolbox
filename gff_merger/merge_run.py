from gff_merger.single_gff_overlap_merger import Single_GFF_Overlap_Merger as gff_mrg
import argparse
import os
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, type=str, help="Input GFF file(s)")
    parser.add_argument("--merge_range", required=True, type=int, help="")
    parser.add_argument("--annotation_type", required=True, type=str, help="")
    args = parser.parse_args()
    input_files = glob.glob(args.gff_in)
    for file in input_files:
        gff_merged = gff_mrg(open(os.path.abspath(file), "r").read(),
                             args.merge_range, args.annotation_type).merge_overlaps()
        print(f"\nWriting output to file: merged_{os.path.basename(file)}")
        outfile = open(f"{os.path.abspath(os.path.join(file, os.pardir))}/merged_{os.path.basename(file)}", "w")
        outfile.write(f"###gff-version 3\n{gff_merged}###")
        outfile.close()
    print("Done!")


main()
