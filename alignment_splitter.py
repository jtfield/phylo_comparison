#! /usr/bin/python3

import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    parser.add_argument('--taxon_num')
    parser.add_argument('-p', action='store_true', help='parsing and splitting of phycorder produced alignment file')
    parser.add_argument('-f', action='store_true', help='parsing and splitting of a genome fasta file')
    return parser.parse_args()

def main():
    args = parse_args()

    seq_file = open(args.align_file, 'r')

    read_file = seq_file.read()

    split_align = read_file.split('>')

    if args.p == True:

        taxon_seq_list = []

        for seq in split_align[1:]:
            split_name_and_seq = seq.split('\n')
            taxon_name = split_name_and_seq[0]
            #print(taxon_name)
            name_split = taxon_name.split('_')
            #print(name_split)
            taxon_number = name_split[1]
            #print(taxon_number)
            if taxon_number == args.taxon_num:
                print("Correct taxon found, beginning sequence split.")

                print(taxon_number)
                n_count = 0
                nuc_count = 0

                seq = split_name_and_seq[1]
                seperate_seq = ''

                n_set = set(['n', 'N'])
                for nuc in seq:
                    if nuc not in n_set:
                        nuc_count+=1
                        seperate_seq+=str(nuc)

                        if nuc_count >= 25:
                            n_count = 0

                    elif nuc in n_set:
                        n_count+=1

                        if n_count >= 50:
                            if nuc_count >= 2:
                                nuc_count = 0
                                taxon_seq_list.append(seperate_seq)
                                seperate_seq = ''

                loci_count = 0
                for small_seq in taxon_seq_list:
                    loci_count+=1
                    loci_file = open(taxon_name + "_" + str(loci_count) + "chunk.fas", "w")
                    loci_file.write('>')
                    loci_file.write(taxon_name)
                    loci_file.write('_' + str(loci_count))
                    loci_file.write('\n')
                    loci_file.write(small_seq)
                    loci_file.write('\n')

                    loci_file.close()

            elif taxon_number != args.taxon_num:
                print("Not the taxon you selected. Skipping.")


        # loci_count = 0
        # for small_seq in taxon_seq_list:
        #     loci_count+=1
        #     loci_file = open(taxon_name + "_" + str(loci_count) + "chunk.fas", "w")
        #     loci_file.write('>')
        #     loci_file.write(taxon_name)
        #     loci_file.write('_' + str(loci_count))
        #     loci_file.write('\n')
        #     loci_file.write(small_seq)
        #     loci_file.write('\n')
        #
        #     loci_file.close()

# for splitting the original genome fasta sequence
    elif args.f == True:
        for seq in split_align[1:]:
            split_name_and_seq = seq.split('\n')
            seq.split('\n', 1)[1]
            # print(seq)

            taxon_name = split_name_and_seq[0]
            # seq = split_name_and_seq[1]
            #
            loci_file = open(taxon_name + "_chunk.fas", "w")
            loci_file.write('>')
            # loci_file.write(taxon_name)
            # loci_file.write('\n')
            loci_file.write(seq)
            loci_file.write('\n')
            #
            loci_file.close()







    seq_file.close()


if __name__ == '__main__':
    main()
# open_file.close()
