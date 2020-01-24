from gff_overlap_merger import GFF_Overlap_Merger as gff_mrg
import argparse
import os
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, type=str, help="Input GFF file(s), comma separated or wildcards")
    parser.add_argument("--gff_out", required=True, type=str, help="Output GFF file.")
    parser.add_argument("--merge_range", required=True, type=int, help="")
    parser.add_argument("--annotation_type", required=True, type=str, help="")
    parser.add_argument("--single_mode", action='store_true', default=False
                        , help="Single file merge mode, if you pass many files, each file will be merged separately")
    args = parser.parse_args()
    input_files = []
    for input_item in args.gff_in.split(','):
        input_files.extend(glob.glob(input_item))
    output_file = os.path.abspath(args.gff_out)
    if args.single_mode:
        for file in input_files:
            gff_merged, count_before, count_after = gff_mrg(open(os.path.abspath(file), "r").read(),
                                                            args.merge_range, args.annotation_type).merge_overlaps()
            print(f"Total {args.annotation_type} after merge: {count_after} of {count_before}")
            print(f"Merged ratio: {round((count_before - count_after) / count_before * 100, 2)}%")
            print(f"\nWriting output to file: {output_file}")
            outfile = open(output_file, "w")
            outfile.write(f"###gff-version 3\n{gff_merged}###")
            outfile.close()
        print("Done!")
    else:
        files_appended = ""
        for file in input_files:
            files_appended += open(os.path.abspath(file), "r").read()
            if not files_appended.endswith('\n'):
                files_appended += "\n"
        files_merged, count_before, count_after = \
            gff_mrg(files_appended, args.merge_range, args.annotation_type).merge_overlaps()
        files_merged.count('\n')
        print(f"Total {args.annotation_type} count after merge: {count_after} of {count_before}")
        print(f"Merged ratio: {round((count_before - count_after) / count_before * 100, 2)}%")
        print(f"\nWriting output to file: {output_file}")
        outfile = open(output_file, "w")
        outfile.write(f"###gff-version 3\n{files_merged}###")
        outfile.close()
        print("Done!")

main()
