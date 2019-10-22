#! /usr/bin/python

import os
import argparse
import re
import json

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast_output_folder')
    #parser.add_argument('--out_file')
    #parser.add_argument('--position_dict_file')
    #parser.add_argument('--suffix')
    return parser.parse_args()


def main():
    args = parse_args()

    duplicate_seperate_dict = {}

    results_list = os.listdir(args.blast_output_folder)
    #print(results_list)

    evalue = "<Hsp_evalue>(.+)</Hsp_evalue>"
    hsp_align_leng = "<Hsp_align-len>.+</Hsp_align-len>"
    gaps = "<Hsp_gaps>.+</Hsp_gaps>"
    query_len = "<Iteration_query-len>.+</Iteration_query-len>"
    hit = "<Hit>"
    file_name_re = "blast_output-(cluster\d+-.fasta-single.fasta)-(cluster\d+-.fasta).txt"

    compile_hit = re.compile(hit)
    compile_evalue = re.compile(evalue)
    compile_hsp_align_leng = re.compile(hsp_align_leng)
    compile_gaps = re.compile(gaps)
    compile_query_len = re.compile(query_len)
    compile_file_name = re.compile(file_name_re)


    
    for file_name in results_list:
        file_open = open(file_name,"r")
        file_read = file_open.read()
        split_count = 0
        hit_findall = re.findall(compile_hit, file_read)
        if hit_findall:
            
            split_file = file_read.split("</Iteration>")
            for split_chunk in split_file:
                split_count+=1
                if split_count == 1 and type(split_chunk) != 'NoneType':
                    #print(split_chunk)
                    evalue_findall = re.findall(compile_evalue, split_chunk)
                    if evalue_findall:
                        if '0' in evalue_findall:
                            check_file_name = re.findall(compile_file_name, file_name)
                            phycord_file = check_file_name[0][0]
                            gon_phy_file = check_file_name[0][1]

                            
                    
                    #hsp_align_leng

        file_open.close()


if __name__ == '__main__':
    main()
