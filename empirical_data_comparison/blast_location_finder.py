#! /usr/bin/env python3
import argparse
import os
import re
import random
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manipulate_seqs_folder')
    parser.add_argument('--long_seqs_folder')
    parser.add_argument('--output_dir')
    # parser.add_argument('--gon_phy', action='store_true', default=False)
    # parser.add_argument('--rapup', action='store_true', default=False)
    return parser.parse_args()

def seq_converter(seq):

    output = {}

    seq_split = seq.split('\n', 1)
    label = seq_split[0]
    seq = seq_split[1]

    len_seq = len(seq_split[1])
    shorter = seq.replace('\n','')

    main_short = Seq(shorter, generic_dna)
    short_comp = main_short.complement()
    short_reverse = main_short[::-1]
    short_rev_comp = main_short.reverse_complement()

    assert len(main_short) > 0
    assert len(short_comp) > 0
    assert len(short_reverse) > 0
    assert len(short_rev_comp) > 0

    assert len(main_short) == len(short_comp) == len(short_reverse) == len(short_rev_comp)
    
    output['original'] = str(main_short)
    output['complement'] = str(short_comp)
    output['reverse'] = str(short_reverse)
    output['reverse_complement'] = str(short_rev_comp)

    return output

def generate_kmer_positions(seq_len):
    output = []
    
    for i in range(0,50):
        kmer_start = random.randint(0,seq_len)
        if(kmer_start + 14 <= seq_len):
            output.append(kmer_start)
    return output


def find_match(long_seq, dict_of_seqs):
    orig_matches = 0
    comp_matches = 0
    rev_matches = 0
    rev_comp_matches = 0
    compare_nums = []
    

    for orientation, seq in dict_of_seqs.items():
        seq_len = len(seq)
        
        get_kmer_starts = generate_kmer_positions(seq_len)

        for position in get_kmer_starts:
            kmer = seq[position:position + 14]
            compile_kmer = re.compile(kmer.upper())
            find_kmer = re.findall(compile_kmer, long_seq.upper())
            if find_kmer:
                if orientation == 'original':
                    orig_matches+=1
                elif orientation == 'complement':
                    comp_matches+=1
                elif orientation == 'reverse':
                    rev_matches+=1
                elif orientation == 'reverse_complement':
                    rev_comp_matches+=1
    
    compare_nums.append(orig_matches)
    compare_nums.append(comp_matches)
    compare_nums.append(rev_matches)
    compare_nums.append(rev_comp_matches)

    print(compare_nums)

    if max(compare_nums) == orig_matches:
        return 'original'
    elif max(compare_nums) == comp_matches:
        return 'complement'
    elif max(compare_nums) == rev_matches:
        return 'reverse'
    elif max(compare_nums) == rev_comp_matches:
        return 'reverse_complement'



def match_long_with_loci(manip_seq_path, long_seq_path, output_dir):
    manip_folder_contents = os.listdir(manip_seq_path)
    long_seqs_folder_contents = os.listdir(long_seq_path)
    file_info_regex = r'-(cluster\d+)--(.+)$'
    file_info_compile = re.compile(file_info_regex)
    long_name_regex = r'single_tax_snippy-core.full.aln-(\w+)$'
    long_name_compile = re.compile(long_name_regex)

    for manip_file in manip_folder_contents:
        # print(manip_file)
        find_info = re.findall(file_info_compile, manip_file)
        if find_info:
            manip_taxon = find_info[0][1]
            manip_locus = find_info[0][0]
            print(manip_taxon)
            print(manip_locus)

            for long_seq in long_seqs_folder_contents:
                # print(long_seq)
                find_long_info = re.findall(long_name_compile, long_seq)
                if find_long_info:
                    long_seq_name = find_long_info[0]
                    if long_seq_name == manip_taxon:
                        # print(manip_taxon)
                        open_manip_file = open(manip_seq_path +'/'+ manip_file,'r')
                        read_manip_file = open_manip_file.read()

                        open_long_seq = open(long_seq_path +'/'+ long_seq, 'r')
                        read_long_seq = open_long_seq.read()

                        convert_manip = seq_converter(read_manip_file)

                        long_seq_split = read_long_seq.split('\n', 1)
                        label = long_seq_split[0]
                        seq = long_seq_split[1]

                        len_seq = len(long_seq_split[1])
                        long_contiguous = seq.replace('\n','')

                        match_maker = find_match(long_contiguous, convert_manip)
                        print(match_maker)

                        output = open(output_dir +'/'+ 'combined-' + manip_locus + '--' + manip_taxon, 'w')
                        output.write(label)
                        output.write('\n')
                        output.write(long_contiguous)
                        output.write('\n')
                        output.write(label)
                        output.write('\n')
                        output.write(convert_manip[match_maker])

                        output.close()
                        open_long_seq.close()
                        open_manip_file.close()





def main():
    args = parse_args()

    path_to_manip_folder = os.path.realpath(args.manipulate_seqs_folder)
    manip_folder_contents = os.listdir(path_to_manip_folder)

    path_to_long_seqs_folder = os.path.realpath(args.long_seqs_folder)
    long_seqs_folder_contents = os.listdir(path_to_long_seqs_folder)

    find_orientation = match_long_with_loci(path_to_manip_folder, path_to_long_seqs_folder, args.output_dir)
    


if __name__ == '__main__':
    main()