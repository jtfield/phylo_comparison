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
        no_hits_found_this_query_count = 0
        for line in xml_file:

            if "<Iteration_query-def>" in line:
                line = line.replace("<Iteration_query-def>","")
                query_name = line.replace("</Iteration_query-def>","")


            elif "<Hit_id>" in line:
                hit = line.replace("<Hit_id>","")
                hit_name = hit.replace("</Hit_id>","")

            elif "<Iteration_message>No hits found</Iteration_message>" in line:
                no_hits_found_this_query_count+=1

            elif "</Iteration>" in line:
                if no_hits_found_this_query_count == 0 :
                    query_name = query_name.strip()
                    hit_name = hit_name.strip()
                    matching_seqs_dict[query_name] = hit_name

                elif no_hits_found_this_query_count >= 1:
                    no_hits_found_this_query_count = 0


            else:
                continue


    file_list = os.listdir(args.dir)
    print(args.dir)
    # print(matching_seqs_dict)

    matched_files_count = 0
    for query_file, match_file in matching_seqs_dict.items():

        query_comp = re.compile("(" + query_file + "chunk.fas)")
        query = re.search(query_comp, str(file_list))

        match_comp = re.compile("(" + match_file + "_chunk.fas)")
        match = re.search(match_comp, str(file_list))

        if query and match:
            query_file_name = (query.group(1))
            match_file_name = (match.group(1))
            matched_files_count+=1

            file_names_list = []
            file_names_list.append(query_file_name)
            file_names_list.append(match_file_name)

            with open(args.dir + "/unaligned_matched_seqs-" + str(matched_files_count) + "-.fa", 'w') as outfile:
                for fname in file_names_list:
                    with open(args.dir + "/" + fname) as infile:
                        for line in infile:
                            outfile.write(line)
            #subprocess.call(['cat ', str(args.dir) + "/" + str(q) , str(args.dir) + "/" + str(m), "> matched_seqs-" + str(matched_files_count) + "-.fa" ],shell=True)

            # combined_files = open("matched_seqs-" + str(matched_files_count) + "-.txt", "w")
            # combined_files.write(query_file_name)
            # combined_files.write(" ")
            # combined_files.write(match_file_name)
            # combined_files.close()
        else:
            print("NO MATCHES")


if __name__ == '__main__':
    main()
