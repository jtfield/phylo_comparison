#! /usr/bin/env python3
import argparse
import os
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder_1')
    parser.add_argument('--folder_2')
    parser.add_argument('--cluster_id_1', default='')
    parser.add_argument('--cluster_id_2', default='')
    parser.add_argument('--output_dir', default='NONE')
    parser.add_argument('--matched_seq_output_dir', default='NONE')
    return parser.parse_args()

def main():
    args = parse_args()

    cluster_id_compile_1 = re.compile(args.cluster_id_1)
    cluster_id_compile_2 = re.compile(args.cluster_id_2)
    
    path_to_input_folder_1 = os.path.realpath(args.folder_1)
    input_folder_contents_1 = os.listdir(path_to_input_folder_1)

    path_to_input_folder_2 = os.path.realpath(args.folder_2)
    input_folder_contents_2 = os.listdir(path_to_input_folder_2)

    matching_folder_1_contents = []
    matching_folder_2_contents = []

    for file_1 in input_folder_contents_1:
        find_cluster = re.findall(cluster_id_compile_1, file_1)
        if find_cluster:
            matching_folder_1_contents.append(file_1)
    
    for file_1 in input_folder_contents_2:
        find_cluster = re.findall(cluster_id_compile_2, file_1)
        if find_cluster:
            matching_folder_2_contents.append(file_1)
        
    print(matching_folder_1_contents)
    print(matching_folder_2_contents)

    cluster_1_taxon_name_regex = ".*-" + args.cluster_id_1 + "--(.+)$"
    compile_cluster_1_name = re.compile(cluster_1_taxon_name_regex)

    cluster_2_taxon_name_regex = ".*-" + args.cluster_id_2 + "--(.+)$"
    compile_cluster_2_name = re.compile(cluster_2_taxon_name_regex)

    for file_1 in matching_folder_1_contents:
        # print(compile_cluster_1_name)
        find_name_1 = re.findall(compile_cluster_1_name, file_1)
        if find_name_1:
            # print(find_name_1)
            for file_2 in matching_folder_2_contents:
                find_name_2 = re.findall(compile_cluster_2_name, file_2)
                if find_name_2:
                    # print(find_name_2)
                    if find_name_1 == find_name_2:
                        print(find_name_1)
                        file_1 = open(path_to_input_folder_1 + '/' + file_1, 'r').read()
                        file_2 = open(path_to_input_folder_2 + '/' + file_2,'r').read()

                        # print(file_1)
                        # print(file_2)
                        if args.output_dir != "NONE":
                            output = open(args.output_dir + "combined_" + args.cluster_id_1 + "_" + args.cluster_id_2 + "_" + ''.join(find_name_1), 'w')
                            output.write(file_1)
                            output.write('\n')
                            output.write(file_2)
                            output.close()

    print("Outputting separate single tax files for second dataset")
    for file_name in matching_folder_2_contents:
        find_name_2 = re.findall(compile_cluster_2_name, file_name)
        if find_name_2:
            print(find_name_2)
            if args.matched_seq_output_dir != "NONE":
                output = open(args.matched_seq_output_dir + '/' + 'single_tax-' + args.cluster_id_2 + '--' + ''.join(find_name_2),'w')
                file_2 = open(path_to_input_folder_2 + '/' + file_name,'r').read()
                output.write(file_2)
                output.close()

if __name__ == '__main__':
    main()