#! /usr/bin/env python3
import argparse
import os
import re
import multiprocessing as mp
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align')
    #parser.add_argument('--align_2')
    parser.add_argument('--output_align_stub', nargs='?', type=str, default="NONE")
    #parser.add_argument('--output_miscalls', nargs='?', type=str, default="NONE")
    #parser.add_argument('--orientation', nargs='?', type=str, default="NONE")
    return parser.parse_args()



def main():
    args = parse_args()

    #pool = mp.Pool(mp.cpu_count())
    print("Number of processors: ", mp.cpu_count())


    align_1 = open(args.align,'r').read()
    #align_2 = open(args.align_2,'r').read()

    align_1_split = align_1.split('\n', 1)
    #align_2_split = align_2.split('\n', 1)

    label = align_1_split[0]
    seq = align_1_split[1]

    len_align_1 = len(align_1_split[1])
    #len_align_2 = len(align_2_split[1])

    
    shorter = seq.replace('\n','')
    #longer = longer.replace('\n','')
    
    #print("Longer sequence length: ", len(longer))
    #print("Shorter sequence length: ", len(shorter))
    

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

   
   
   #Produce reverse sequence
    output_file = open(args.output_align_stub + "-reverse.fasta", 'w+')
    output_file.write(label)
    output_file.write('\n')
    output_file.write(str(short_reverse))
    output_file.close()


   #Produce compliment sequence
    output_file = open(args.output_align_stub + "-complement.fasta", 'w+')
    output_file.write(label)
    output_file.write('\n')
    output_file.write(str(short_comp))
    output_file.close()


    #Produce reverse complement sequence
    output_file = open(args.output_align_stub + "-reverse_complement.fasta", 'w+')
    output_file.write(label)
    output_file.write('\n')
    output_file.write(str(short_rev_comp))
    output_file.close()



    #print("Number of processors: ", mp.cpu_count())

if __name__ == '__main__':
    main()
