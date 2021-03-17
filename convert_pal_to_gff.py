import os
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pal_in", required=True, help="", type=str)
    args = parser.parse_args()
    with open(os.path.abspath(args.pal_in), "r") as f:
        f_str = f.read()
        f_str_rec_l = f_str.split("Palindromes of:  ")
        f_str_rec_l = f_str_rec_l[1:]
        all_pos = []
        for rec in f_str_rec_l:
            rec_split = rec.split("Palindromes:")
            seq_id = rec_split[0].split(" \n", maxsplit=1)[0]
            pos = get_pos(rec_split[1])
            all_pos.extend([make_gff_record(x, seq_id) for x in pos])
        with open(os.path.abspath(f"{args.pal_in}.gff"), "w") as of:
            of.write("\n".join(all_pos))


def get_pos(str_in):
    r_lst = []
    str_in_lst = str_in.split("\n")
    counter = 0
    for indx, line in enumerate(str_in_lst):
        if "|" in line:
            counter += 1
            line1_lst = list(filter(None, str_in_lst[indx - 1].split(" ")))
            line2_lst = list(filter(None, str_in_lst[indx + 1].split(" ")))
            pos_dict = {'id': f"{counter}",
                        'stem_loop_start': line1_lst[0],
                        'stem_loop_end': line2_lst[0],
                        'stem_loop_len': int(line2_lst[0]) - int(line1_lst[0]) + 1,
                        'stem_len': int(line1_lst[2]) - int(line1_lst[0]) + 1,
                        'loop_len': int(line2_lst[2]) - int(line1_lst[2]) - 1}
            r_lst.append(pos_dict)
    return r_lst


def make_gff_record(in_dict, seqid):
    return f"{seqid}\tEMBOSS_Palindrome\tinverted_repeat\t{in_dict['stem_loop_start']}\t{in_dict['stem_loop_end']}\t.\t+\t.\t"\
           f"ID=IR_{seqid}_{in_dict['id']}_f;Name=Inverted_Repeat_{seqid}_{in_dict['id']}_f;"\
           f"stem_loop_len={in_dict['stem_loop_len']};stem_len={in_dict['stem_len']};loop_len={in_dict['loop_len']}\n"\
           f"{seqid}\tEMBOSS_Palindrome\tinverted_repeat\t{in_dict['stem_loop_start']}\t{in_dict['stem_loop_end']}\t.\t-\t.\t"\
           f"ID=IR_{seqid}_{in_dict['id']}_r;Name=Inverted_Repeat_{seqid}_{in_dict['id']}_r;"\
           f"stem_loop_len={in_dict['stem_loop_len']};stem_len={in_dict['stem_len']};loop_len={in_dict['loop_len']}"\

main()