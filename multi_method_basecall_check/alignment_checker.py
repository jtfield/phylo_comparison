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



def main():
    args = parse_args()

    align_1 = open(args.align_1,'r').read()
    align_2 = open(args.align_2,'r').read()

    len_align_1 = len(align_1)
    len_align_2 = len(align_2)

    combine_list = [len_align_1, len_align_2]

    longer = ''
    shorter = ''
    if len_align_1 == max(combine_list):
        longer = align_1
        shorter = align_2

    elif len_aling_2 == max(combine_list):
        longer = align_2
        shorter = align_1

    # make sure the right files are assigned to the longer and shorter variables
    assert len(longer) > len(shorter)

    indentical_bases_by_position = {}
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
            print(combined_positions)




if __name__ == '__main__':
    main()
