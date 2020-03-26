from pattern_locations_fetcher import PatternLocationsFetcher as plf
from wiggle_merger import WiggleMerger as wm
from Bio import SeqIO
import glob
import argparse


def main():
    # Params
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_in", required=True, help="RefSeq fasta file", type=str)
    parser.add_argument("--wigs_in", required=True,
                        help="Coverage wiggle file(s), Must contain forward and reverse files", type=str)
    parser.add_argument("--gff_out", required=True, help="GFF output file name for terminators", type=str)
    parser.add_argument("--pre_signal_offset", required=True, help="", type=int)
    parser.add_argument("--post_signal_offset", required=True, help="", type=int)
    parser.add_argument("--merge_range", default=30, required=False, help="", type=int)
    parser.add_argument("--max_interruption", default=3, required=False, help="", type=int)
    parser.add_argument("--peak_percentage", required=True, type=int,
                        help="")

    parser.add_argument("--base", default="T", required=False, help="", type=str)
    parser.add_argument("--min_base_count", default=4, required=False, help="", type=int)
    args = parser.parse_args()
    fasta_parsed = SeqIO.parse(glob.glob(args.fasta_in)[0], "fasta")
    wig_files = glob.glob(args.wigs_in)
    f_wigs_max, r_wigs_max = wm(wig_files, fasta_parsed).merge_wigs_by_max()
    ret_list = []
    counters = {}
    # The following line is repeated due to the previous iterator exhaustion
    fasta_parsed = SeqIO.parse(glob.glob(args.fasta_in)[0], "fasta")
    x = 1
    for seq_record in fasta_parsed:
        f_seq_str = str(seq_record.seq)
        accession = seq_record.id
        f_positions = plf(f_seq_str, "f", args.base, args.max_interruption).fetch_locations()
        #r_positions = plf(f_seq_str, "r", args.base, args.max_interruption).fetch_locations()
        f_positions = [x for x in f_positions if len(x) >= args.min_base_count]
        for index, pos in enumerate(f_positions):
            peak_info = get_peaks(f_wigs_max[accession],
                                  [pos[0] - args.pre_signal_offset,
                                   pos[-1] + args.post_signal_offset])
            if peak_info is None:
                counters[f'pos_not_in_cov_{accession}_f'] = 0
                continue
            if not pos[0] <= peak_info['peak_pos'] <= pos[-1]:
                if pos[-1] < peak_info['peak_pos']:
                    pos[-1] = peak_info['peak_pos']
                if pos[0] > peak_info['peak_pos']:
                    pos[0] = peak_info['peak_pos']
            f_positions[index] = [pos[0], pos[-1],
                                  peak_info['peak_pos'],
                                  peak_info['peak_score'],
                                  peak_info['peak_percentage']]
        for i in f_positions:
            print(i)



def get_peaks(wig_df, pos):
    wig_df = wig_df[wig_df[0].between(pos[0], pos[-1])]
    if wig_df.empty:
        return None
    else:
        peak_index = wig_df[1].idxmax()
        peak_pos = wig_df.at[peak_index, 0]
        peak_score = wig_df.at[peak_index, 1]
        #wig_df.duplicated(keep='last')
        wig_df = wig_df.drop([peak_index])
        peak_percentage = round((peak_score - wig_df[1].mean()) / peak_score * 100, 2)
        return {'peak_pos': peak_pos,
                'peak_score': peak_score,
                'peak_percentage': peak_percentage}

def merge_interval_lists(list_in, merge_range):
    merge_range += 2
    list_out = []
    for loc in list_in:
        if len(list_out) == 0:
            list_out.append(loc)
        else:
            if loc[0] in range(list_out[-1][0], list_out[-1][-1] + merge_range):
                list_out[-1][-1] = loc[-1]
            else:
                list_out.append(loc)
    return list_out


main()
