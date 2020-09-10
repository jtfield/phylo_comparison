#! /usr/bin/python
# program to split a multiple sequence alignment into seperate sequence files for each taxon

import os
import argparse
import re
import csv

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_file_dir')
    #parser.add_argument('--position_csv_file')
    parser.add_argument('--concatenated_fasta')
    return parser.parse_args()


def main():
    args = parse_args()

    input_file = open(args.concatenated_fasta, 'r')
    read_input = input_file.read()
    len_read_input = len(read_input)
    #print(len_read_input)

    # make sure this isnt an empty file
    assert(len_read_input >= 1)

    # split file into multiple sequences
    split_taxa = read_input.split('>')
    #print(split_taxa)

    
    file_name = args.concatenated_fasta
    file_name = os.path.basename(file_name)
    file_name = file_name.replace('.fasta','')
    file_name = file_name.replace('.fas','')
    file_name = file_name.replace('.aln', '')

    for chunk in split_taxa:
        if len(chunk) < 1:
            #print(chunk)
            print("empty chunk")
        else:
            name_split = chunk.split('\n', 1)
            #print(name_split)
            #print(len(name_split))
            # test that you got two things from this split: a name and a sequence
            assert(len(name_split) == 2)
            name = name_split[0]
            name = name.replace('.fasta','')
            #name = name.replace('.fas','')
            name = name.replace('.fq', '')
            name = name.replace('.fastq','')
            name = name.replace('.fas','')
            seq = name_split[1]
            new_file_name = file_name + '_' + name + '.fas'
            print(new_file_name)

            # make sure file name is a real sort of file name
            assert(len(new_file_name) > 1)

            # open file, add sequence and carry on
            output = args.out_file_dir + '/' + new_file_name
            open_output = open(output, 'w')
            open_output.write('>')
            open_output.write(chunk)
            open_output.write('\n')
            open_output.close()




if __name__ == '__main__':
    main()

