#! /usr/bin/env python3
import argparse
import os
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder')
    # parser.add_argument('--cluster_id', default='')
    parser.add_argument('--output_file')
    return parser.parse_args()


def main():
    args = parse_args()

    cluster_id_compile = re.compile("single_tax_gon_phy-cluster[0-9]+")

    path_to_input_folder = os.path.realpath(args.input_folder)
    input_folder_contents = os.listdir(path_to_input_folder)

    print(input_folder_contents)

    best_blast_scores_per_cluster = {}
    cluster_name_set = []
    for file_name in input_folder_contents:
        find_cluster_id = re.findall(cluster_id_compile, file_name)
        if find_cluster_id:
            print(find_cluster_id)
            print(type(find_cluster_id))
            print(type(find_cluster_id[0]))

            cluster_name_set.append(find_cluster_id[0])

    cluster_name_set = list(set(cluster_name_set))

    # print(cluster_name_set)

    find_alignment = '<Hsp_midline>(.+)</Hsp_midline>'
    compile_find_alignment = re.compile(find_alignment)

    best_matches = {}
    for cluster_id in cluster_name_set:
        blast_scores_for_cluster = {}
        cluster_id = cluster_id + '-'
        compile_id = re.compile(cluster_id)
        print(cluster_id)
        for file in input_folder_contents:
            # print(file)
            find_id = re.findall(compile_id, file)
            if find_id:
                first_hit = 0
                with open(path_to_input_folder + '/' + file) as data:
                    for line in data:
                        align_finder = re.findall(compile_find_alignment, line)
                        if align_finder:
                            if first_hit == 0:
                                # print(line)

                                for alignment in align_finder:
                                    total_length = 0
                                    correct_align = 0
                                    for nuc_line in alignment:
                                        if nuc_line == '|':
                                            correct_align+=1
                                        total_length+=1
                                    blast_scores_for_cluster[file] = correct_align

                            first_hit+=1
        best_match_for_locus = max(blast_scores_for_cluster, key=blast_scores_for_cluster.get)
        best_matches[cluster_id] = best_match_for_locus

    print(best_matches)

    best_match_locus = '<BlastOutput_db>(.+)</BlastOutput_db>'
    compile_locus = re.compile(best_match_locus)
    
    locus_match_dict = {}
    for key, value in best_matches.items():
        with open(path_to_input_folder + '/' + value) as best_match_file:
            for line in best_match_file:
                find_locus_file = re.findall(compile_locus, line)
                if find_locus_file:
                    locus_match_dict[key] = ''.join(find_locus_file)

    print(locus_match_dict)

    # output = open(args.output_file, 'w')
    # output.write()
    locus_count = 0
    for key, value in locus_match_dict.items():
        locus_count+=1
        output = open(args.output_file + str(locus_count) + '.txt', 'w')
        output.write(key)
        output.write('\n')
        output.write(value)
        # output.write('\n')

        output.close() 
        
        
    print("DONE")

if __name__ == '__main__':
    main()