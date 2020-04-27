#! /usr/bin/env python3
import argparse
import os
import re
import multiprocessing as mp
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_1')
    parser.add_argument('--align_2')
    parser.add_argument('--align_3')
    parser.add_argument('--align_4')
    #parser.add_argument('--output_align_stub', nargs='?', type=str, default="NONE")
    #parser.add_argument('--output_miscalls', nargs='?', type=str, default="NONE")
    #parser.add_argument('--orientation', nargs='?', type=str, default="NONE")
    return parser.parse_args()

def alignment_fixer(read_file):
    split_seqs = read_file.split('>')

    output = []
    seq_1 = []
    seq_2 = []
    seq_count = 0
    for seq_and_name in split_seqs:
        if len(seq_and_name) > 0:
            seq_count+=1
            split_seq_and_name = seq_and_name.split('\n', 1)
            #print(split_seq_and_name)
            name = split_seq_and_name[0]
            seq = split_seq_and_name[1]
        
            fixed_seq = seq.replace('\n','')
            if seq_count == 1:
                seq_1.append(name)
                seq_1.append(fixed_seq)
            elif seq_count == 2:
                seq_2.append(name)
                seq_2.append(fixed_seq)

    output.append(seq_1)
    output.append(seq_2)

    return output

def nuc_counter(sequence):
    nuc_count = 0
    for nuc in sequence[1]:
        if nuc != '-':
            nuc_count+=1
    
    return nuc_count

def trim_gaps(short_seq):
    output = []
    leading_gaps = 0
    trailing_gaps = 0
    
    leading_gap_match = "^(-*)\w"
    trailing_gap_match = "\w(-*)$"
    compile_leading_match = re.compile(leading_gap_match)
    compile_trailing_match = re.compile(trailing_gap_match)
    find_leading = re.findall(compile_leading_match, short_seq)
    find_trailing = re.findall(compile_trailing_match, short_seq)

    if compile_leading_match:
        #print("found leading")
        #print(find_leading)
        leading_gaps = len(find_leading[0])
        print(leading_gaps)
    
    if compile_trailing_match:
        #print("found trailing")
        #print(find_trailing)
        trailing_gaps = len(find_trailing[0])
        print(leading_gaps)

    output.append(leading_gaps)
    output.append(trailing_gaps)
    
    return output
    

def comparison(list_of_list_of_seqs):

    seq_1 = None
    seq_2 = None
    seq_count = 0
    for list_of_seqs in list_of_list_of_seqs:
        seq_count+=1
        if seq_count == 1:
            seq_1 = list_of_seqs
        elif seq_count == 2:
            seq_2 = list_of_seqs

    #print(seq_1)
    #print(seq_2)

    count_nucs_1 = nuc_counter(seq_1)
    count_nucs_2 = nuc_counter(seq_2)
    
    #print(count_nucs_1)
    #print(count_nucs_2)
    
    nuc_lens = [count_nucs_1, count_nucs_2]
    small_seq = min(nuc_lens)

    shorter = None
    longer = None
    if small_seq == count_nucs_1:
        shorter = seq_1
        longer = seq_2

    elif small_seq == count_nucs_2:
        shorter = seq_2
        longer = seq_1

    #print(shorter)
    #print(longer)
    
    shorter_seq = shorter[1]
    longer_seq = longer[1]
    
    get_gaps = trim_gaps(shorter_seq) 
    #print(get_gaps)

    trimmed_shorter = shorter_seq[get_gaps[0]:-get_gaps[1]]
    #print(trimmed_shorter)

    trimmed_longer = longer_seq[get_gaps[0]:-get_gaps[1]]
    

    

def main():
    args = parse_args()

    #pool = mp.Pool(mp.cpu_count())
    print("Number of processors: ", mp.cpu_count())


    align_1 = open(args.align_1,'r').read()
    #align_2 = open(args.align_2,'r').read()
    #align_3 = open(args.align_3,'r').read()
    #align_3 = open(args.align_4,'r').read()

    parse_align_1 = alignment_fixer(align_1)
    #print(parse_align_1)

    compare_seqs = comparison(parse_align_1) 
    print(compare_seqs)

    #align_1_split = align_1.split('\n', 1)
    #align_2_split = align_2.split('\n', 1)
    #align_3_split = align_3.split('\n', 1)
    #align_4_split = align_4.split('\n', 1)

    #label = align_1_split[0]
    #seq = align_1_split[1]

    #len_align_1 = len(align_1_split[1])
    #len_align_2 = len(align_2_split[1])

    
    #shorter = seq.replace('\n','')
    #longer = longer.replace('\n','')
    
    #print("Longer sequence length: ", len(longer))
    #print("Shorter sequence length: ", len(shorter))
    
    #print(len(shorter))



    #main_short = Seq(shorter, generic_dna)
    #short_comp = main_short.complement()
    #short_reverse = main_short[::-1]
    #short_rev_comp = main_short.reverse_complement()
   
   
   #Produce reverse sequence
    #output_file = open(args.output_align_stub + "-reverse.fasta", 'w')
    #output_file.write(label)
    #output_file.write('\n')
    #output_file.write(str(short_reverse))
    #output_file.close()


   #Produce compliment sequence
    #output_file = open(args.output_align_stub + "-complement.fasta", 'w')
    #output_file.write(label)
    #output_file.write('\n')
    #output_file.write(str(short_comp))
    #output_file.close()


    #Produce reverse complement sequence
    #output_file = open(args.output_align_stub + "-reverse_complement.fasta", 'w')
    #output_file.write(label)
    #output_file.write('\n')
    #output_file.write(str(short_rev_comp))
    #output_file.close()



    #print("Number of processors: ", mp.cpu_count())

if __name__ == '__main__':
    main()
