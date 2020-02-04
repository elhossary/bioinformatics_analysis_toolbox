from gff_overlap_merger import GFF_Overlap_Merger as gff_mrg
import argparse
import os
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, type=str, nargs='+',
                        help="Input GFF file(s) --gff_in file 1, --gff_in file 2")
    parser.add_argument("--gff_out", required=False, type=str, help="Output GFF file.")
    parser.add_argument("--merge_range", required=False, type=int, default=0,
                        help="Distance between 2 sequences to merge in one if they are not overlapping")
    parser.add_argument("--annotation_type", required=False, type=str, help="", default="Annotation type after merge")
    parser.add_argument("--single_mode", action='store_true', default=False,
                        help="Single file merge mode, if you pass many files, each file will be merged separately")
    parser.add_argument("--annotate", default='all', required=False, choices=['all', 'overlaps', 'no_overlaps'],
                        help="choose one: all, overlaps, or no_overlaps")
    args = parser.parse_args()
    if not args.single_mode and not args.gff_out:
        parser.error("ERROR:  --gff_out argument is required when --single_mode argument id not used.")
    input_files = []
    for item in args.gff_in:
        input_files.extend(glob.glob(item))
    input_files = list(set(input_files))
    if args.single_mode:
        args.annotation_type = ""
        for file in input_files:
            if args.annotate == 'overlaps':
                prefix = "overlaps_"
            elif args.annotate == 'no_overlaps':
                prefix = "no_overlaps_"
            else:
                prefix = "merged_"
            output_base_name = f"{prefix}{os.path.basename(file)}"
            output_path = os.path.abspath(os.path.join(file, os.pardir))
            output_file = f"{output_path}/{output_base_name}"

            gff_merged, count_before, count_after = gff_mrg(open(os.path.abspath(file), "r").read(),
                                                            args.annotation_type, args.merge_range,
                                                            args.annotate).merge_overlaps()
            print(f"Total annotations count before merge:\t{count_before}")
            print(f"Total annotations count after merge:\t{count_after}")
            print(f"Difference:\t{count_before - count_after}")
            if args.annotate == 'overlaps':
                overlap_ratio = round(count_after / count_before * 100, 2)
            elif args.annotate == 'no_overlaps':
                overlap_ratio = \
                    round((count_before - (count_after + (count_before - count_after) / 2)) / count_before * 100, 2)
            else:
                overlap_ratio = round((count_before - count_after) / count_before * 100, 2)
            print(f"Overlap ratio: {overlap_ratio}%")
            print(f"Writing output to file: {output_file}")
            outfile = open(output_file, "w")
            outfile.write(f"###gff-version 3\n{gff_merged}###")
            outfile.close()
        print("Done!\n")
    else:
        output_file = os.path.abspath(args.gff_out)
        files_appended = ""
        for file in input_files:
            files_appended += open(os.path.abspath(file), "r").read()
            if not files_appended.endswith('\n'):
                files_appended += "\n"
        files_merged, count_before, count_after = \
            gff_mrg(files_appended, args.annotation_type, args.merge_range, args.annotate).merge_overlaps()
        files_merged.count('\n')
        print(f"Total {args.annotation_type} count before merge:\t{count_before}")
        print(f"Total {args.annotation_type} count after merge:\t{count_after}")
        print(f"Difference:\t{count_before - count_after}")
        if args.annotate == 'overlaps':
            overlap_ratio = round(count_after / count_before * 100, 2)
        elif args.annotate == 'no_overlaps':
            overlap_ratio = \
                round((count_before - (count_after + (count_before - count_after) / 2)) / count_before * 100, 2)
        else:
            overlap_ratio = round((count_before - count_after) / count_before * 100, 2)

        print(f"Decrease ratio: {overlap_ratio}%")
        print(f"Writing output to file: {output_file}")
        outfile = open(output_file, "w")
        outfile.write(f"###gff-version 3\n{files_merged}###")
        outfile.close()
        print("Done!\n")

main()