#! /usr/bin/env python3
import argparse
import os
import re
import csv
import numpy
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--loci_dir')
    parser.add_argument('--file_suffix')
    return parser.parse_args()

def loci_len_check(path_to_loci, loci_folder_contents, suffix):

    loci_lengths = []
    for num, locus_file in enumerate(loci_folder_contents):
        if locus_file.endswith(suffix):
            file_lines = 0
            input_file = open(path_to_loci + '/' + locus_file)
            read_input = input_file.read()
            split_input = read_input.split('>')
            for taxon in split_input:
                split_name_and_seq = taxon.split('\n', 1)
                if len(split_name_and_seq) > 1:
                    file_lines+=1
                    if file_lines == 1:
                        loci_lengths.append(len(split_name_and_seq[1]))
                    elif file_lines == 2:
                        break


            input_file.close()
                
        
    
    return loci_lengths

def organize_and_describe(list_of_lengths):
    output_values = {}
    loci_lens_dict = {}
    min_len = min(list_of_lengths)
    max_len = max(list_of_lengths)
    num_loci = len(list_of_lengths)
    avg_len = sum(list_of_lengths) / len(list_of_lengths)
    output_values['min'] = min_len
    output_values['max'] = max_len
    output_values['avg'] = avg_len
    output_values['num'] = num_loci
    # output_values['lengths'] = list_of_lengths
    print(output_values)

    loci_lens_dict['lengths'] = list_of_lengths

    df = pd.DataFrame(loci_lens_dict)

    df.to_csv('loci_lengths.csv')



def main():
    args = parse_args()

    path_to_align_folder = os.path.realpath(args.loci_dir)
    align_folder_contents = os.listdir(path_to_align_folder)

    process_loci = loci_len_check(path_to_align_folder, align_folder_contents, args.file_suffix)

    report_values = organize_and_describe(process_loci)



if __name__ == '__main__':
    main()