#! /usr/bin/python3

import sys
import os
import argparse
import re
import subprocess

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq_set_1_dir')
    parser.add_argument('--seq_set_2_dir')
    parser.add_argument('--ref_seq_dir')
    parser.add_argument('--ref_suffix')
    return parser.parse_args()

def main():
    args = parse_args()

    seq_set_1_list = os.listdir(args.seq_set_1_dir)
    seq_set_2_list = os.listdir(args.seq_set_2_dir)
    ref_seq_dir_list = os.listdir(args.ref_seq_dir)

    #print(ref_seq_dir_list)

    # make sure the directories actually have files in them
    assert(len(seq_set_1_list) > 0)
    assert(len(seq_set_2_list) > 0)
    assert(len(ref_seq_dir_list) > 0)

    tax_reg = '(taxon_\d+)'
    cluster_reg = '(cluster\d+-)'

    compile_tax_reg = re.compile(tax_reg)
    compile_cluster_reg = re.compile(cluster_reg)

    #print(seq_set_2_list)
    # test that file names have been formatted properly
    seq_1_join = ''.join(seq_set_1_list)
    seq_2_join = ''.join(seq_set_2_list)
    ref_join = ''.join(ref_seq_dir_list)
    


    seq_1_check = re.search(compile_tax_reg, seq_1_join)
    seq_2_check = re.search(compile_tax_reg, seq_2_join)
    seq_ref_check = re.search(compile_tax_reg, ref_join)
    assert(type(seq_1_check.group()) == str)
    assert(type(seq_2_check.group()) == str)
    assert(type(seq_ref_check.group()) == str)


    taxon_count = 0
    # Begin matching sequence files to reference files
    for file_name in seq_set_1_list:
        if file_name.endswith('.fas'):
            file_match_list = []
            taxon_count+=1
            seq_1_search = re.search(compile_tax_reg, file_name)
            cluster_check = re.search(compile_cluster_reg, file_name)
            cluster_name = cluster_check.group()
            if seq_1_search:
                #print(seq_1_search.group())
                #cluster_check = re.search(compile_cluster_reg, file_name)
                #cluster_name = cluster_check.group()
                current_taxon = seq_1_search.group() + '_'
                #print(current_taxon)
                #print(type(current_taxon))
                compile_current_taxon = re.compile(current_taxon)
                for ref_name in ref_seq_dir_list:
                    if ref_name.endswith(args.ref_suffix):
                        ref_match = re.search(compile_current_taxon, ref_name)
                        if ref_match:
                            #print(ref_name)
                            file_match_list.append(args.seq_set_1_dir + '/' + file_name)
                            file_match_list.append(args.ref_seq_dir + '/' + ref_name)
            
            # make sure you only have 2 file names in the list before writing to file
            assert(len(file_match_list) == 2)

            # write matched file names to files so they can be blasted against eachother
            output = open(args.seq_set_1_dir + '/' + cluster_name + current_taxon + '_' + 'matched_files_.txt', 'w')
            output.write(file_match_list[0])
            output.write('\n')
            output.write(file_match_list[1])
            output.close()
    


    for file_name in seq_set_2_list:
        if file_name.endswith('.fas'):
            file_match_list = []
            taxon_count+=1
            seq_2_search = re.search(compile_tax_reg, file_name)
            cluster_check = re.search(compile_cluster_reg, file_name)
            cluster_name = cluster_check.group()
            if seq_2_search:
                #print(seq_1_search.group())
                current_taxon = seq_2_search.group() + '_'
                #print(current_taxon)
                #print(type(current_taxon))
                compile_current_taxon = re.compile(current_taxon)
                for ref_name in ref_seq_dir_list:
                    if ref_name.endswith(args.ref_suffix):
                        ref_match = re.search(compile_current_taxon, ref_name)
                        if ref_match:
                            #print(ref_name)
                            file_match_list.append(args.seq_set_2_dir + '/' + file_name)
                            file_match_list.append(args.ref_seq_dir + '/' + ref_name)

            # make sure you only have 2 file names in the list before writing to file
            assert(len(file_match_list) == 2)

            # write matched file names to files so they can be blasted against eachother
            output = open(args.seq_set_2_dir + '/' + cluster_name + current_taxon + '_' + 'matched_files_.txt', 'w')
            output.write(file_match_list[0])
            output.write('\n')
            output.write(file_match_list[1])
    print('finished file matching')




if __name__ == '__main__':
    main()

