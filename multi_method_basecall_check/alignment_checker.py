#! /usr/bin/env python3
import argparse
import os
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_1')
    parser.add_argument('--align_2')
    return parser.parse_args()


def split_seq(string):
    split_string = list(string)

    return split_string

def check_alignment(list_of_paired_nucs):
    identical_nucs = 0
    for pair in list_of_paired_nucs:
        
        assert len(pair) == 2
        if pair[0] == pair[1]:
            identical_nucs+=1

    return identical_nucs


def main():
    args = parse_args()

    align_1 = open(args.align_1,'r').read()
    align_2 = open(args.align_2,'r').read()

    align_1_split = align_1.split('\n', 1)
    align_2_split = align_2.split('\n', 1)

    #print(align_1_split)
    #print(align_2_split)

    len_align_1 = len(align_1_split[1])
    len_align_2 = len(align_2_split[1])

    combine_list = [len_align_1, len_align_2]

    longer = ''
    shorter = ''
    longer_label = ''
    shorter_label = ''
    if len_align_1 == max(combine_list):
        longer = align_1_split[1]
        longer_label = align_1_split[0]
        shorter = align_2_split[1]
        shorter_label = align_2_split[0]

    elif len_aling_2 == max(combine_list):
        longer = align_2_split[1]
        longer_label = align_2_split[0]
        shorter = align_1_split[1]
        shorter_label = align_1_split[0]

    # make sure the right files are assigned to the longer and shorter variables
    assert len(longer) > len(shorter)
    
    print(len(longer))
    print(len(shorter))

    indentical_bases_by_position = 0
    identical_bases = 0
    split_short = split_seq(shorter)
    for num in range(0, len(longer)):
        #longer_segment = longer[num:len(shorter)]
        longer_segment = longer[num:num + len(shorter)]
        if len(longer_segment) >= len(shorter):
        
            split_long = split_seq(longer_segment)
            #print(len(split_long))
            #print(len(longer))
            
            # make sure that the split subsequence of the longer sequence is not shorter than the shorter sequence (ugh)
            assert len(split_long) >= len(shorter)
           
            #zip both lists together by elements
            combined_positions = list(map(list, zip(split_long, split_short)))
            #print(combined_positions)
            check_identity = check_alignment(combined_positions)
            #print(check_identity)

            if check_identity > identical_bases:
                identical_bases = check_identity
                identical_bases_by_position = num

    print(identical_bases)
    print(identical_bases_by_position)

    output_file = open('test_output_align.fasta','w')
    output_file.write(longer_label)
    output_file.write('\n')
    output_file.write(longer)
    output_file.write('\n')
    output_file.write(shorter_label)
    output_file.write('\n')
    output_file.write((identical_bases_by_position * '-') + shorter)
    output_file.close()

if __name__ == '__main__':
    main()
