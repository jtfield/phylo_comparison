#! /usr/bin/env python3
import argparse
import os
import re
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder_1')
    parser.add_argument('--folder_2')
    parser.add_argument('--cluster_id_1', default='')
    parser.add_argument('--cluster_id_2', default='')
    parser.add_argument('--output_dir', default='NONE')
    #parser.add_argument('--matched_seq_output_dir', default='NONE')
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
    taxa_names = []

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

    assert len(matching_folder_1_contents) > 0
    assert len(matching_folder_2_contents) > 0

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
                        taxa_names.append(find_name_1[0])

                        file_1 = open(path_to_input_folder_1 + '/' + file_1, 'r').read()
                        file_2 = open(path_to_input_folder_2 + '/' + file_2,'r').read()

                        align_2_split = file_2.split('\n', 1)

                        label = align_2_split[0]
                        seq = align_2_split[1]

                        len_align_2 = len(align_2_split[1])
                        shorter = seq.replace('\n','')

                        main_short = Seq(shorter)
                        short_comp = main_short.complement()
                        short_reverse = main_short[::-1]
                        short_rev_comp = main_short.reverse_complement()

                        main_short = str(main_short)
                        short_comp = str(short_comp)
                        short_reverse = str(short_reverse)
                        short_rev_comp = str(short_rev_comp)

                        assert len(main_short) > 0
                        assert len(short_comp) > 0
                        assert len(short_reverse) > 0
                        assert len(short_rev_comp) > 0

                        print("SEQ_LENGTHS")
                        print(len(main_short))
                        print(len(short_comp))
                        print(len(short_reverse))
                        print(len(short_rev_comp))

                        # print(file_1)
                        # print(file_2)
                        if args.output_dir != "NONE":
                            output = open(args.output_dir + "combined_original_" + args.cluster_id_1 + "_" + args.cluster_id_2 + "--" + ''.join(find_name_1) + '--', 'w')
                            output.write(file_1)
                            output.write('\n')
                            output.write(label)
                            output.write('\n')
                            output.write(main_short)
                            output.close()

                            output = open(args.output_dir + "combined_reverse_" + args.cluster_id_1 + "_" + args.cluster_id_2 + "--" + ''.join(find_name_1) + '--', 'w')
                            output.write(file_1)
                            output.write('\n')
                            output.write(label)
                            output.write('\n')
                            output.write(short_reverse)
                            output.close()

                            output = open(args.output_dir + "combined_complement_" + args.cluster_id_1 + "_" + args.cluster_id_2 + "--" + ''.join(find_name_1) + '--', 'w')
                            output.write(file_1)
                            output.write('\n')
                            output.write(label)
                            output.write('\n')
                            output.write(short_comp)
                            output.close()

                            output = open(args.output_dir + "combined_reverse_complement_" + args.cluster_id_1 + "_" + args.cluster_id_2 + "--" + ''.join(find_name_1) + '--', 'w')
                            output.write(file_1)
                            output.write('\n')
                            output.write(label)
                            output.write('\n')
                            output.write(short_rev_comp)
                            output.close()

    # taxa_names = list(set(taxa_names))
    assert len(taxa_names) > 1
    names_output = open(args.output_dir + 'taxa_name_list.txt','w')
    for taxon_name in list(set(taxa_names)):
        names_output.write(taxon_name)
        names_output.write('\n')
    names_output.close()

    #print("Outputting separate single tax files for second dataset")
    #for file_name in matching_folder_2_contents:
    #    find_name_2 = re.findall(compile_cluster_2_name, file_name)
    #    if find_name_2:
    #        print(find_name_2)
    #        if args.matched_seq_output_dir != "NONE":
    #            output = open(args.matched_seq_output_dir + '/' + 'single_tax-' + args.cluster_id_2 + '--' + ''.join(find_name_2),'w')
    #            file_2 = open(path_to_input_folder_2 + '/' + file_name,'r').read()
    #            output.write(file_2)
    #            output.close()

if __name__ == '__main__':
    main()
