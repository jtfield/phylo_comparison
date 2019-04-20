#! /usr/bin/python3

import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    parser.add_argument('--taxon_num')
    return parser.parse_args()

def main():
    args = parse_args()

    seq_file = open(args.align_file, 'r')

    read_file = seq_file.read()

    split_align = read_file.split('>')

    for seq in split_align[1:]:
        split_name_and_seq = seq.split('\n')
        taxon_name = split_name_and_seq[0]
        name_split = taxon_name.split('_')
        taxon_number = name_split[1]
        if taxon_number == args.taxon_num:
            print(taxon_number)
            seq = split_name_and_seq[1]

            # make a complete alignment

            # Split each seq into seperate parts
            split_seq = seq.split('n'*100)

            # start counting each part of the sequence and open files to get those sequences
            loci_count = 0
            for small_seq in split_seq:
                loci_count+=1
                loci_file = open(taxon_name + "_loci_" + str(loci_count) + ".fas", "w")
                loci_file.write('>')
                loci_file.write(taxon_name)
                loci_file.write('_' + str(loci_count))
                loci_file.write('\n')
                loci_file.write(small_seq)
                loci_file.write('\n')

                loci_file.close()









    seq_file.close()


if __name__ == '__main__':
    main()
# open_file.close()
