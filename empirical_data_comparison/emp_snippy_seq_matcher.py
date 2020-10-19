#! /usr/bin/env python3
import argparse
import os
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manipulate_seqs_folder')
    parser.add_argument('--ref_seqs_folder')
    parser.add_argument('--output_align_stub')
    # parser.add_argument('--cluster_id_1', default='')
    # parser.add_argument('--cluster_id_2', default='')
    parser.add_argument('--align_output_dir')
    parser.add_argument('--list_output_dir')
    # parser.add_argument('--matched_seq_output_dir')
    return parser.parse_args()

def main():
    args = parse_args()

    cluster_and_name_regex = '-(cluster\d+)--(.+)$'
    compile_regex = re.compile(cluster_and_name_regex)

    path_to_input_folder_1 = os.path.realpath(args.manipulate_seqs_folder)
    input_folder_contents_1 = os.listdir(path_to_input_folder_1)

    path_to_input_folder_2 = os.path.realpath(args.ref_seqs_folder)
    input_folder_contents_2 = os.listdir(path_to_input_folder_2)

    print(input_folder_contents_1)
    print(input_folder_contents_2)
    print(path_to_input_folder_1)
    print(path_to_input_folder_2)

    list_of_loci = []
    list_of_taxa = []

    for file_name in input_folder_contents_1:
        find_info = re.findall(compile_regex, file_name)
        if find_info:
            print(find_info)
            assert len(find_info) == 1
            assert len(find_info[0]) == 2
            cluster = find_info[0][0]
            tax_name = find_info[0][1]

            align = open(path_to_input_folder_1 + '/' + file_name,'r')
            align_1 = align.read()

            align_1_split = align_1.split('\n', 1)

            label = align_1_split[0]
            seq = align_1_split[1]

            len_align_1 = len(align_1_split[1])
            shorter = seq.replace('\n','')

            main_short = Seq(shorter, generic_dna)
            short_comp = main_short.complement()
            short_reverse = main_short[::-1]
            short_rev_comp = main_short.reverse_complement()

            assert len(main_short) > 0
            assert len(short_comp) > 0
            assert len(short_reverse) > 0
            assert len(short_rev_comp) > 0

            print("SEQ_LENGTHS")
            print(len(main_short))
            print(len(short_comp))
            print(len(short_reverse))
            print(len(short_rev_comp))

            # FIND THE MATCHING TAXON SEQUENCE FILE IN THE SECOND DIRECTORY
            compile_tax_name = re.compile('-' + tax_name + '$')
            for file_name_2 in input_folder_contents_2:
                # print(file_name_2)
                find_name = re.findall(compile_tax_name, file_name_2)
                if find_name:
                    
                    list_of_loci.append(cluster)
                    list_of_taxa.append(tax_name)

                    matching_file_2 = open(path_to_input_folder_2 + '/' + file_name_2,'r')
                    # file_1_read = matching_file_1.read()
                    file_2_read = matching_file_2.read()

                    #Produce normal sequence
                    output_file = open(args.align_output_dir + '/' + args.output_align_stub + '_'+ cluster + '--' + tax_name + "--original.fasta", 'w+')
                    output_file.write(label)
                    output_file.write('\n')
                    output_file.write(str(main_short))
                    output_file.write('\n')
                    output_file.write(file_2_read)
                    output_file.close()

                    #Produce reverse sequence
                    output_file = open(args.align_output_dir + '/' + args.output_align_stub + '_'+ cluster + '--' + tax_name + "--reverse.fasta", 'w+')
                    output_file.write(label)
                    output_file.write('\n')
                    output_file.write(str(short_reverse))
                    output_file.write('\n')
                    output_file.write(file_2_read)
                    output_file.close()


                    # Produce compliment sequence
                    output_file = open(args.align_output_dir + '/' + args.output_align_stub + '_'+ cluster + '--' + tax_name + "--complement.fasta", 'w+')
                    output_file.write(label)
                    output_file.write('\n')
                    output_file.write(str(short_comp))
                    output_file.write('\n')
                    output_file.write(file_2_read)
                    output_file.close()


                    #Produce reverse complement sequence
                    output_file = open(args.align_output_dir + '/' + args.output_align_stub + '_'+ cluster + '--' + tax_name + "--reverse_complement.fasta", 'w+')
                    output_file.write(label)
                    output_file.write('\n')
                    output_file.write(str(short_rev_comp))
                    output_file.write('\n')
                    output_file.write(file_2_read)
                    output_file.close()

    output_names = open(args.list_output_dir + '/' + "taxa_list.txt", 'w')
    for name in list(set(list_of_taxa)):
        output_names.write(name)
        output_names.write('\n')
    output_names.close

    output_taxa = open(args.list_output_dir + '/' + "loci_list.txt", 'w')
    for locus in list(set(list_of_loci)):
        output_taxa.write(locus)
        output_taxa.write('\n')
    output_taxa.close

    # print(input_folder_contents_1)
    # print(input_folder_contents_2)
    # print(path_to_input_folder_1)
    # print(path_to_input_folder_2)
    # for file_name in input_folder_contents_1:
    #     find_info = re.findall(compile_regex, file_name)
    #     if find_info:
    #         print(find_info)
    #         assert len(find_info) == 1
    #         assert len(find_info[0]) == 2
    #         cluster = find_info[0][0]
    #         tax_name = find_info[0][1]
    #         compile_tax_name = re.compile(tax_name + '$')
    #         for file_name_2 in input_folder_contents_2:
    #             # print(file_name_2)
    #             find_name = re.findall(compile_tax_name, file_name_2)
    #             if find_name:
    #                 # file_1_path = path_to_input_folder_1 + '/' + file_name
    #                 # file_2_path = path_to_input_folder_2 + '/' + file_name_2
    #                 # print(file_1_path)
    #                 # print(file_2_path)

    #                 matching_file_1 = open(path_to_input_folder_1 + '/' + file_name,'r')
    #                 matching_file_2 = open(path_to_input_folder_2 + '/' + file_name_2,'r')
    #                 file_1_read = matching_file_1.read()
    #                 file_2_read = matching_file_2.read()
    #                 output = open(args.output_dir + '/' + 'combined_seqs_' + cluster + '_' + tax_name,'w')
    #                 output.write(file_1_read)
    #                 output.write(file_2_read)
    #                 output.close()
    #                 matching_file_1.close()
    #                 matching_file_2.close()






if __name__ == '__main__':
    main()