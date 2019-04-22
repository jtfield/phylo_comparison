#! /usr/bin/python3

import sys
import os
import argparse
import re

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


    print(matching_seqs_dict)



if __name__ == '__main__':
    main()
