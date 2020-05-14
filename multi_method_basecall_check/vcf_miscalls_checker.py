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
    for name in main_vcf_header:
        print(name)
        name_compile = re.compile('(' + name.strip() + ')')
        name_match = re.search(name_compile, simulated_vcf_file)
        if name_match:
            print("FOUND")
            print(name_match)


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
    #for line in combined_lines_frame:
    #    print(line)
    #    print(len(line))

    df = pd.DataFrame(combined_lines_frame, columns=get_head)
    #print(df)
    
    print(args.vcf_dir)
    find_correct_name = find_taxon(get_head, args.vcf_dir)




if __name__ == '__main__':
    main()
