import argparse
import pandas as pd
from os import path
from Bio import SeqIO
from more_itertools import consecutive_groups


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff_in", required=True, help="", type=str)
    parser.add_argument("--refseq_in", required=True, help="", type=str)
    parser.add_argument("--tss_in", required=False, help="", type=str)
    parser.add_argument("--energy_values_in", required=True, help="", type=str)
    parser.add_argument("--end_range", default=10, required=False, help="", type=int)
    parser.add_argument("--offset", default=0, required=False, help="", type=int)
    parser.add_argument("--gff_out", required=True, help="", type=str)
    args = parser.parse_args()
    fasta_parsed = SeqIO.parse(path.abspath(args.refseq_in), "fasta")
    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(path.abspath(args.gff_in), names=col_names, sep="\t", comment="#")
    str_out = ""

    f_seq = ""
    r_seq = ""
    step_factor = 0
    step_height = 0
    t_count = 0
    longest_t_count = 0
    energy_values_list = parse_energy_values_file(args.energy_values_in)
    if args.tss_in is not None:
        tss_df = pd.read_csv(path.abspath(args.tss_in), names=col_names, sep="\t", comment="#")
    else:
        tss_df = None
    for seq_record in fasta_parsed:
        for index, row in gff_df.iterrows():
            if row['seqid'] == seq_record.id:
                f_seq = str(seq_record.seq)
                r_seq = str(seq_record.reverse_complement().seq)[::-1]  # reverse of reverse_complement
            else:
                continue
            energy_value = ""
            if row['strand'] == "+":
                seq = f_seq[int(row['end']) - args.end_range - 1:int(row['end']) + args.offset]
                t_count = seq.count("T")
                longest_t_count = get_longest_continuous_t(seq)
                for item in energy_values_list:
                    if row['seqid'] in item[0] and "+" in item[0] and \
                            parse_attributes(row['attributes'])['id'] in item[0] and \
                            seq.replace("T", "U") in item[1]:
                        energy_value = item[2].split(' ')[-1].replace('(', '').replace(')', '')
                        break
            elif row['strand'] == "-":
                seq = r_seq[int(row['start']) - args.offset - 1:int(row['start']) + args.end_range]
                t_count = seq.count("T")
                longest_t_count = get_longest_continuous_t(seq)
                for item in energy_values_list:
                    if row['seqid'] in item[0] and "-" in item[0] and \
                            parse_attributes(row['attributes'])['id'] in item[0] and \
                            seq.replace("T", "U")[::-1] in item[1]:
                        energy_value = item[2].split(' ')[-1].replace('(', '').replace(')', '')
                        break
            else:
                print("Fatal error")
            if tss_df is not None:
                tmp = tss_df[(tss_df["seqid"] == row['seqid']) &
                             (tss_df["strand"] == row['strand']) &
                             (tss_df["start"].between(row['start'], row['end']))]["attributes"].values.tolist()

                if len(tmp) == 1:
                    step_factor = parse_attributes(tmp[0])["ave_step_factor"]
                    step_height = parse_attributes(tmp[0])["ave_step_height"]
                elif len(tmp) == 0:
                    step_factor = 0
                    step_height = 0
                else:
                    step_factor, step_height = get_strongest_tss_from_attrib(tmp)

            str_out += \
                f"{row['seqid']}\t" + \
                f"{row['source']}\t" + \
                f"{row['type']}\t" + \
                f"{int(row['start'])}\t" + \
                f"{int(row['end'])}\t" + \
                f"{row['score']}\t" + \
                f"{row['strand']}\t" + \
                f"{row['phase']}\t" + \
                f"{row['attributes']}" + \
                f";seq_len={int(row['end']) - int(row['start']) + 1}" + \
                f";T_count_in_last_{args.end_range}_bases_with_{args.offset}_offset={t_count}" + \
                f";longest_continuous_T={longest_t_count}" + \
                f";energy_value={energy_value}" + \
                f";ave_step_height={step_height}" + \
                f";ave_step_factor={step_factor}" + \
                "\n"
    print("Writing GFF file...")
    outfile = open(path.abspath(args.gff_out), "w")
    outfile.write(f"{str_out}")
    outfile.close()


def get_longest_continuous_t(sequence_str):
    count = 0
    indices = [i for i, a in enumerate(sequence_str, 1) if a == "T"]
    indices.sort()
    signals = [list(group) for group in consecutive_groups(indices)]
    for i in signals:
        if len(i) > count:
            count = len(i)
    return count


def parse_energy_values_file(path_str):
    energy_values_file = open(path.abspath(path_str), "r")
    counter = 0
    energy_values_list = []
    tmp = ""
    for line in energy_values_file.readlines():
        counter += 1
        tmp += line
        if counter == 3:
            energy_values_list.append(tmp.split("\n"))
            tmp = ""
            counter = 0
    return energy_values_list


def get_strongest_tss_from_attrib(lst_in):
    lst_in = [parse_attributes(item) for item in lst_in]
    return max([float(x['ave_step_factor']) for x in lst_in]), max([float(x['ave_step_height']) for x in lst_in])


def parse_attributes(attr_str):
    return {k.lower(): v for k, v in dict(item.split("=") for item in attr_str.split(";")).items()}

main()