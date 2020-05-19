#! /usr/bin/env python3
import argparse
import os
import re
import numpy as np
import pandas as pd
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--main_vcf')
    parser.add_argument('--vcf_dir')
    parser.add_argument('--msa')
    return parser.parse_args()

def get_header(file_name):
    header = []
    with open(file_name) as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                if not line.startswith('##'):
                    header.append(line.strip().split())
    header = header[0]
    return header

def find_taxon(main_vcf_header, simulated_vcf_file):
    taxon_name = ''
    for name in main_vcf_header:
        #print(name)
        name_compile = re.compile('(' + name.strip() + ')')
        name_match = re.search(name_compile, simulated_vcf_file)
        if name_match:
            #print("FOUND")
            #print(name_match)
            taxon_name = name

    return taxon_name




def check_miscalls(prime_dataframe, secondary_dataframe):
    output_values = []
    for index, row in prime_dataframe.iterrows():
        position = row['POS']
        ref = row['REF']
        alt = row['ALT']
        #print(secondary_dataframe['POS'] == position)

def main():
    args = parse_args()

    get_head = get_header(args.main_vcf)
    #print(get_head)

    #head_df = pd.Dataframe
    lines_frame = []
    with open(args.main_vcf) as main_vcf:
        line_count_1 = 0
        line_count_2 = 0
        for num, line in enumerate(main_vcf):
            if not line.startswith('#'):
                line_count_1+=1
                line = line.split()
                #print(get_head)
                #print(len(get_head))
                #print(line)
                #print(len(line))
                #print("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ")
                lines_frame.append(line)
    #print(lines_frame)

    combined_lines_frame = []
    #print(len(lines_frame))
    temp_line = ''
    for num in range(len(lines_frame)):
        if num % 2 == 0:
            temp_line = lines_frame[num]
        else:
            #print(num)
            #print(lines_frame[num])
            combined_lines_frame.append(temp_line + lines_frame[num])

    #print(combined_lines_frame)
    vcf_path = args.vcf_dir + "/snps.vcf"
    find_correct_name = find_taxon(get_head, args.vcf_dir)
    taxon_name_compile = re.compile(find_correct_name)

    print(find_correct_name)
    print(len(find_correct_name))

    #print(get_head)
    pos_position = 0
    ref_position = 0
    alt_position = 0
    tax_of_interest_position = 0
    for num, info in enumerate(get_head):
        #print(info)
        #print(type(info))
        info = info.strip('\n')
        if info == 'POS':
            pos_position = num
        elif info == 'REF':
            ref_position = num
        elif info == 'ALT':
            alt_position = num
        elif info == find_correct_name:
            tax_of_interest_position = num

    #print(pos_position)
    #print(ref_position)
    #print(alt_position)
    #print(tax_of_interest_position)

    prime_vcf_dict_tax_of_interest = {}
    for line in combined_lines_frame:
        assert len(line) == len(get_head)
        prime_vcf_dict_tax_of_interest[line[pos_position]] = [line[ref_position], line[alt_position], line[tax_of_interest_position]]

    #print(prime_vcf_dict_tax_of_interest)

    get_single_tax_head = get_header(vcf_path)

    single_tax_lines = []
    with open(vcf_path) as single_tax_vcf:
        for num, line in enumerate(single_tax_vcf):
            if not line.startswith('#'):
                single_tax_lines.append(line)


    sing_tax_pos_position = 0
    sing_tax_ref_position = 0
    sing_tax_alt_position = 0
    #tax_of_interest_position = 0
    for num, info in enumerate(get_single_tax_head):
        #print(info)
        #print(type(info))
        info = info.strip('\n')
        if info == 'POS':
            sing_tax_pos_position = num
        elif info == 'REF':
            sing_tax_ref_position = num
        elif info == 'ALT':
            sing_tax_alt_position = num
    
    #print(get_single_tax_head)
    #print(len(get_single_tax_head))
    #print(single_tax_lines[0])

    sing_tax_vcf_dict_tax_of_interest = {}
    for line in single_tax_lines:
        line = line.split()
        #print(line)
        #print(len(line))
        assert len(line) == len(get_single_tax_head)
        sing_tax_vcf_dict_tax_of_interest[line[sing_tax_pos_position]] = [line[sing_tax_ref_position], line[sing_tax_alt_position]] 

    print(sing_tax_vcf_dict_tax_of_interest)
    #print(prime_vcf_dict_tax_of_interest)

    should_have_variants = 0
    actually_has_variants = 0
    has_variant_dict = {}
    should_have_variant_dict = {}

    for key, value in prime_vcf_dict_tax_of_interest.items():
        #print(key)
        if value[2] != '0':
            should_have_variants+=1
            should_have_variant_dict[key] = value
            if key in sing_tax_vcf_dict_tax_of_interest:
                actually_has_variants+=1
                has_variant_dict[key] = value

    #print(has_variant_dict)
    #print("waffle")
    #print(should_have_variant_dict)


    for key, value in has_variant_dict.items():
        should_have_variant_dict.pop(key, None)

    #print(should_have_variants)
    print(should_have_variant_dict)

    #print(actually_has_variants)
    print(has_variant_dict)

    
    msa = open(args.msa, 'r')
    read_msa = msa.read()

    split_msa = read_msa.split('>')

    n_positions = []
    gap_positions = []
    sequence_of_interest = []
    for taxon in split_msa:
        if len(taxon) > 0:
            split_taxon = taxon.split('\n',1)
            taxon_name_search = re.search(taxon_name_compile, split_taxon[0])
            if taxon_name_search:
                #print(split_taxon[0])
                contig_seq = ''.join(split_taxon[1].split('\n'))
                #print(contig_seq)
                for num, letter in enumerate(contig_seq):
                    sequence_of_interest.append(letter)
                    if letter == 'N':
                        n_positions.append(num)
                    elif letter == '-':
                        gap_positions.append(num)

#    #print(n_positions)
#    #print(gap_positions)
#
#    non_missing_data_miscalls = {}
#    for key, value in doesnt_have_variant_dict.items():
#        key = int(key)
#        if key not in n_positions:
#            if key not in gap_positions:
#                non_missing_data_miscalls[key] = value
#
#    #print(non_missing_data_miscalls)
#
#    for key, value in non_missing_data_miscalls.items():
#        #print(key)
#        #print(value)
#        #print(sequence_of_interest[key-1])
#        continue


    #df = pd.DataFrame(combined_lines_frame, columns=get_head)
    
    #find_correct_name = find_taxon(get_head, args.vcf_dir)
    #print(find_correct_name)

    #vcf_path = args.vcf_dir + "/snps.vcf"
    #get_single_tax_head = get_header(vcf_path)
    
    #print(get_single_tax_head)

    #single_tax_df = pd.read_csv(vcf_path, comment = '#', sep='\t', names=get_single_tax_head)
    #print(single_tax_df)
    
    #single_tax_values_of_interest = single_tax_df[['POS', 'REF', 'ALT']].copy()
    #print(single_tax_values_of_interest)
    
    #specific_tax_in_original_vcf = df[['POS', 'REF', 'ALT' ,find_correct_name]].copy()

    #orig_to_dict = specific_tax_in_original_vcf.to_dict()
    #print(orig_to_dict)

    #single_tax_to_dict = single_tax_values_of_interest.to_dict()
    #print(single_tax_to_dict)


    #check = check_miscalls(specific_tax_in_original_vcf, single_tax_values_of_interest)

if __name__ == '__main__':
    main()
