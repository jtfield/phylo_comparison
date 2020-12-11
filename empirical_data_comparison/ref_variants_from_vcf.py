#! /usr/bin/env python3
import argparse
import os
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    parser.add_argument('--output_align')
    return parser.parse_args()

def parse_vcf(path_to_vcf):
    output_nucleotides = []
    output_name = []
    output = []

    with open(path_to_vcf, 'r') as vcf:
        vcf_lines = vcf.readlines()
        for line in vcf_lines:
            if not line.startswith("#"):
                split_line = line.split("\t")
                if len(output_name) < 1:
                    output_name.append(split_line[0])
                if split_line[3] != '.':
                    output_nucleotides.append(split_line[3])
    output_nucleotides = ''.join(output_nucleotides)
    output.append(output_name)
    output.append(output_nucleotides)
    return output
                



def main():
    args = parse_args()

    vcf = os.path.realpath(args.vcf_file)

    parser = parse_vcf(vcf)

    output = open(args.output_align, 'w')
    output.write('>' + parser[0][0])
    output.write('\n')
    output.write(parser[1])

if __name__ == '__main__':
    main()