#! /usr/bin/python3

import sys
import os
import argparse
import re
import subprocess

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast_xml_file')
    parser.add_argument('--dir')
    return parser.parse_args()

def main():
    args = parse_args()

    # establish dictionary for matching sequences together
    matching_seqs_dict = {}
    query_name = ""
    hit_name = ""
    with open(args.blast_xml_file, 'r') as xml_file:
        loci_count = 0
        for line in xml_file:

            if "<Iteration_query-def>" in line:
                line = line.replace("<Iteration_query-def>","")
                query_name = line.replace("</Iteration_query-def>","")


            elif "<Hit_id>" in line:
                hit = line.replace("<Hit_id>","")
                hit_name = hit.replace("</Hit_id>","")

            elif "</Iteration>" in line:
                query_name = query_name.strip()
                hit_name = hit_name.strip()
                matching_seqs_dict[query_name] = hit_name

            else:
                continue


    file_list = os.listdir(args.dir)

    # print(matching_seqs_dict)

    matched_files_count = 0
    for query_file, match_file in matching_seqs_dict.items():

        query_comp = re.compile("(" + query_file + "chunk.fas)")
        query = re.search(query_comp, str(file_list))

        match_comp = re.compile("(" + match_file + "_chunk.fas)")
        match = re.search(match_comp, str(file_list))

        if query and match:
            q = (query.group(1))
            m = (match.group(1))
            matched_files_count+=1
            subprocess.call(['cat ', str(q) , str(q), "> matched_seqs-" + str(matched_files_count) + "-.fa" ],shell=True)




if __name__ == '__main__':
    main()
