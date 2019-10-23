#! /usr/bin/python

import os
import argparse
import re
import json

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa_folder')
    parser.add_argument('--out_file')
    parser.add_argument('--position_dict_file')
    #parser.add_argument('--suffix')
    return parser.parse_args()


def main():
    args = parse_args()

    msa_list = os.listdir(args.msa_folder)
    print(msa_list)

    json_data = open(args.position_dict_file,'r')
    read_dict = json_data.read()
    pos_dict = json.loads(read_dict)
    print(pos_dict)

    loci_for_use = []
    loci_count = 0
    for pos, locus_name_and_len in pos_dict.items():
        for leng, name in locus_name_and_len.items():
            #print(type(leng))
            if int(leng) >= 2000 and loci_count < 10 and name in msa_list:
                loci_for_use.append(str(name))
                loci_count+=1
     
    print(loci_for_use)

    open_output = open(args.out_file,"w")
    for name in loci_for_use:
        open_output.write(name)
        open_output.write("\n")
    open_output.close()

if __name__ == '__main__':
    main()