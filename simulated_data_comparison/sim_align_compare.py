#! /usr/bin/env python3
import argparse
import os
import re
import multiprocessing as mp
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_1')
    parser.add_argument('--align_2')
    parser.add_argument('--align_3')
    parser.add_argument('--align_4')
    parser.add_argument('--output_stub', nargs='?', type=str, default="NONE")
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

    #output.append(identical_nucs)
    #output.append(non_identical_nucs)
    #output.append(gaps)
    #output.append(identical_positions)
    #output.append(non_identical_positions)
    #output.append(gap_positions) 
    output.append(leading_gaps)
    output.append(trailing_gaps)

    return output
    
def check_alignment(list_of_paired_nucs):
    total_nucs = 0
    identical_nucs = 0
    non_identical_nucs = 0
    gaps = 0
    gap_positions = []
    identical_positions = []
    non_identical_positions = []
    gap_set = ['-','N']
    nuc_set = ['A','C','G','T']
    output = []
    for num, pair in enumerate(list_of_paired_nucs):
        assert len(pair) == 2
        total_nucs+=1
        #print(pair[0])
        #print(pair[1])
        if str(pair[0].upper()) in gap_set or str(pair[1].upper()) in gap_set:
            #if pair[0] in gap_set:
            #    print(pair[0])
            #elif pair[1] in gap_set:
            #    print(pair[1])
            gaps+=1
            gap_positions.append(num)
            #print('gap')

        elif pair[0] and pair[1] not in gap_set:
            #print('no gap')
            if pair[0].upper() == pair[1].upper():
                identical_nucs+=1
                identical_positions.append(num)

            elif pair[0].upper() != pair[1].upper():
                non_identical_nucs+=1
                non_identical_positions.append(num)
   
    output.append(identical_nucs)
    output.append(non_identical_nucs)
    output.append(gaps)
    output.append(total_nucs)
    output.append(identical_positions)
    output.append(non_identical_positions)
    output.append(gap_positions)


#    identical_nucs = 0
#    non_identical_nucs = 0
#    nuc_set = []
#    for pair in list_of_paired_nucs:
#
#        assert len(pair) == 2
#        if pair[0].upper() == pair[1].upper():
#            identical_nucs+=1
#        elif pair[0].upper() != pair[1].upper():
#            non_identical_nucs+=1
#    nuc_set.append(identical_nucs)
#    nuc_set.append(non_identical_nucs)

    #return identical_nucs
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
   
    #print(len(shorter_seq))
    #print(len(longer_seq))

    len_short = nuc_counter(shorter)
    len_long = nuc_counter(longer)

    #print(shorter_seq)
    #print(longer_seq)

    get_gaps = trim_gaps(shorter_seq) 
    #print(get_gaps)

    trimmed_shorter = ''
    trimmed_longer = ''
    print("waffle1")
    if get_gaps[1] != 0:

        #print(shorter_seq[get_gaps:-get_gaps[1]])
        trimmed_shorter = shorter_seq[get_gaps[0]:-get_gaps[1]]
        trimmed_longer = longer_seq[get_gaps[0]:-get_gaps[1]]
        #print(trimmed_shorter)
        #print(trimmed_longer)
    elif get_gaps[1] == 0:
        #print(shorter_seq[get_gaps[0]:])
        trimmed_shorter = shorter_seq[get_gaps[0]:]
        trimmed_longer = longer_seq[get_gaps[0]:]
        #print(trimmed_shorter)
        #print(trimmed_longer)
    print("waffle2")

    #trimmed_longer = longer_seq[get_gaps[0]:-get_gaps[1]]


    #trimmed_shorter = shorter_seq[get_gaps[0]:-get_gaps[1]]
    #print(trimmed_shorter)

    #trimmed_longer = longer_seq[get_gaps[0]:-get_gaps[1]]
    
    split_trimmed_short = list(trimmed_shorter)
    split_trimmed_long = list(trimmed_longer)

    combined_positions = list(map(list, zip(split_trimmed_long, split_trimmed_short)))
    analyze_alignment = check_alignment(combined_positions)

    #print(analyze_alignment)
    analyze_alignment.append(len_short)
    analyze_alignment.append(len_long)
    
    return analyze_alignment

def main():
    args = parse_args()

    #pool = mp.Pool(mp.cpu_count())
    print("Number of processors: ", mp.cpu_count())


    align_1 = open(args.align_1,'r').read()
    align_2 = open(args.align_2,'r').read()
    align_3 = open(args.align_3,'r').read()
    align_4 = open(args.align_4,'r').read()

    parse_align_1 = alignment_fixer(align_1)
    #print(parse_align_1)
    parse_align_2 = alignment_fixer(align_2)
    parse_align_3 = alignment_fixer(align_3)
    parse_align_4 = alignment_fixer(align_4)

    compare_seqs_1 = comparison(parse_align_1)
    compare_seqs_2 = comparison(parse_align_2)
    compare_seqs_3 = comparison(parse_align_3)
    compare_seqs_4 = comparison(parse_align_4)
    print("align 1", compare_seqs_1[0])
    print("align 1", compare_seqs_1[1])
    print("align 1", compare_seqs_1[2])
    print("align 1", compare_seqs_1[3])
    print("align 2", compare_seqs_2[0])
    print("align 2", compare_seqs_2[1])
    print("align 2", compare_seqs_2[2])
    print("align 2", compare_seqs_2[3])
    print("align 3", compare_seqs_3[0])
    print("align 3", compare_seqs_3[1])
    print("align 3", compare_seqs_3[2])
    print("align 3", compare_seqs_3[3])
    print("align 4", compare_seqs_4[0])
    print("align 4", compare_seqs_4[1])
    print("align 4", compare_seqs_4[2])
    print("align 4", compare_seqs_4[3])

    compare_identical_nucs = [compare_seqs_1[0], compare_seqs_2[0], compare_seqs_3[0], compare_seqs_4[0]]
    best_align = max(compare_identical_nucs)

    output_file = open(args.output_stub, 'w')
    
    if best_align == compare_seqs_1[0]:
        output_file.write(args.align_1)
        output_file.write("\n")
        output_file.write(">identical_nucs")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[0]))
        output_file.write("\n")
        output_file.write(">non_identical_nucs")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[1]))
        output_file.write("\n")
        output_file.write(">gaps")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[2]))
        output_file.write("\n")
        output_file.write(">total_nucleotides")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[3]))
        output_file.write("\n")
        output_file.write(">identical_nucs_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[4]))
        output_file.write("\n")
        output_file.write(">non_identical_nucs_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[5]))
        output_file.write("\n")
        output_file.write(">gaps_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[6]))
        output_file.write("\n")
        output_file.write(">unadjusted_short_seq_length")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[7]))
        output_file.write("\n")
        output_file.write(">unadjusted_long_seq_length")
        output_file.write("\n")
        output_file.write(str(compare_seqs_1[8]))


    elif best_align == compare_seqs_2[0]:
        output_file.write(args.align_2)
        output_file.write("\n")
        output_file.write(">identical_nucs")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[0]))
        output_file.write("\n")
        output_file.write(">non_identical_nucs")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[1]))
        output_file.write("\n")
        output_file.write(">gaps")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[2]))
        output_file.write("\n")
        output_file.write(">total_nucleotides")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[3]))
        output_file.write("\n")
        output_file.write(">identical_nucs_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[4]))
        output_file.write("\n")
        output_file.write(">non_identical_nucs_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[5]))
        output_file.write("\n")
        output_file.write(">gaps_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[6]))
        output_file.write("\n")
        output_file.write(">unadjusted_short_seq_length")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[7]))
        output_file.write("\n")
        output_file.write(">unadjusted_long_seq_length")
        output_file.write("\n")
        output_file.write(str(compare_seqs_2[8]))

    elif best_align == compare_seqs_3[0]:
        output_file.write(args.align_3)
        output_file.write("\n")
        output_file.write(">identical_nucs")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[0]))
        output_file.write("\n")
        output_file.write(">non_identical_nucs")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[1]))
        output_file.write("\n")
        output_file.write(">gaps")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[2]))
        output_file.write("\n")
        output_file.write(">total_nucleotides")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[3]))
        output_file.write("\n")
        output_file.write(">identical_nucs_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[4]))
        output_file.write("\n")
        output_file.write(">non_identical_nucs_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[5]))
        output_file.write("\n")
        output_file.write(">gaps_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[6]))
        output_file.write("\n")
        output_file.write(">unadjusted_short_seq_length")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[7]))
        output_file.write("\n")
        output_file.write(">unadjusted_long_seq_length")
        output_file.write("\n")
        output_file.write(str(compare_seqs_3[8]))

    elif best_align == compare_seqs_4[0]:
        output_file.write(args.align_4)
        output_file.write("\n")
        output_file.write(">identical_nucs")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[0]))
        output_file.write("\n")
        output_file.write(">non_identical_nucs")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[1]))
        output_file.write("\n")
        output_file.write(">gaps")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[2]))
        output_file.write("\n")
        output_file.write(">total_nucleotides")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[3]))
        output_file.write("\n")
        output_file.write(">identical_nucs_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[4]))
        output_file.write("\n")
        output_file.write(">non_identical_nucs_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[5]))
        output_file.write("\n")
        output_file.write(">gaps_positions")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[6]))
        output_file.write("\n")
        output_file.write(">unadjusted_short_seq_length")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[7]))
        output_file.write("\n")
        output_file.write(">unadjusted_long_seq_length")
        output_file.write("\n")
        output_file.write(str(compare_seqs_4[8]))


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
