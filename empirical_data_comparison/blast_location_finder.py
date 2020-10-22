#! /usr/bin/env python3
import argparse
import os
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manipulate_seqs_folder')
    parser.add_argument('--long_seqs_folder')
    parser.add_argument('--blast_files_folder')
    parser.add_argument('--cluster_id')
    parser.add_argument('--output_dir')
    # parser.add_argument('--gon_phy', action='store_true', default=False)
    # parser.add_argument('--rapup', action='store_true', default=False)
    return parser.parse_args()

def check_alignment_direction(blast_files, blast_folder_path, cluster):
    cluster_dict = {}
    align_numbers = []
    hit_regex = r'<Hit>'
    hit_compile = re.compile(hit_regex)

    start_regex = r'<Hsp_hit-from>(\d+)</Hsp_hit-from>'
    start_compile = re.compile(start_regex)

    end_regex = r'<Hsp_hit-to>(\d+)</Hsp_hit-to>'
    end_compile = re.compile(end_regex)

    cluster_regex = cluster + '-'
    compile_regex = re.compile(cluster_regex)
    
    for file_name in blast_files:
        
        find_cluster = re.findall(compile_regex, file_name)
        if find_cluster:
            current_file = open(blast_folder_path +'/'+ file_name, 'r')
            read_file = current_file.read()
            find_hit = re.findall(hit_compile, read_file)
            if find_hit:
                file_start_stop = []
                find_start = re.findall(start_compile, read_file)
                find_stop = re.findall(end_compile, read_file)
                if find_start:
                    file_start_stop.append(find_start[0])
                if find_stop:
                    file_start_stop.append(find_stop[0])
                # print(file_start_stop)
                assert len(file_start_stop) == 2
                align_numbers.append(file_start_stop)
    cluster_dict[cluster] = align_numbers
    
    # print(cluster_dict)
    return cluster_dict      

def analyze_direction(dict_of_align_nums):
    strandedness = ''
    increasing = 0
    decreasing = 0
    for key, value in dict_of_align_nums.items():
        for start_stop in value:
            start = int(start_stop[0])
            stop = int(start_stop[1])
            if start > stop:
                decreasing+=1
            elif start < stop:
                increasing+=1
    if increasing > 0 and increasing > decreasing:
        strandedness = 'left_to_right'
    elif decreasing > 0 and increasing < decreasing:
        strandedness = 'right_to_left'

    return strandedness    


def find_long_seq_size(align_nums, direction):
    seq_size = []
    starts = []
    stops = []
    start = 0
    stop = 0
    for key, value in align_nums.items():
        for start_stop in value:
            potential_start = int(start_stop[0])
            potential_stop = int(start_stop[1])
            starts.append(potential_start)
            stops.append(potential_stop)
    
    if direction == 'left_to_right':
        start = min(starts)
        stop = max(stops)
    elif direction == 'right_to_left':
        start = max(starts)
        stop = min(stops)
    seq_size.append(start)
    seq_size.append(stop)
   
    return seq_size

def convert_short_seq(direction, file_contents):
    split_file = file_contents.split('\n',1)
    assert len(split_file) > 1
    name = split_file[0]
    seq = split_file[1]

    if direction == 'right_to_left':
        reverse = seq[::-1]
        output = name + '\n' + reverse
        return output
    elif direction == 'left_to_right':
        return file_contents

def trim_long_seq(start, stop, file_contents):
    split_file = file_contents.split('\n',1)
    assert len(split_file) > 1
    name = split_file[0]
    seq = split_file[1]

    seq = ''.join(seq)

    padded_start = start - 100
    padded_stop = stop + 100

    seq = seq[padded_start:padded_stop]

    output = name + '\n' + seq
    return output

def resize_long(sizes, direction, long_seqs, manip_seqs, outdir):
    start = min(sizes)
    stop = max(sizes)
    manip_folder_contents = os.listdir(manip_seqs)
    long_seqs_folder_contents = os.listdir(long_seqs)
    file_info_regex = r'-(cluster\d+)--(.+)$'
    file_info_compile = re.compile(file_info_regex)
    long_name_regex = r'single_tax_snippy-core.full.aln-(\w+)$'
    long_name_compile = re.compile(long_name_regex)

    for manip_file in manip_folder_contents:
        # print(manip_file)
        find_info = re.findall(file_info_compile, manip_file)
        if find_info:
            manip_taxon = find_info[0][1]
            manip_locus = find_info[0][0]
            for long_seq in long_seqs_folder_contents:
                # print(long_seq)
                find_long_info = re.findall(long_name_compile, long_seq)
                if find_long_info:
                    long_seq_name = find_long_info[0]
                    if long_seq_name == manip_taxon:
                        # print(manip_taxon)
                        open_manip_file = open(manip_seqs +'/'+ manip_file,'r')
                        read_manip_file = open_manip_file.read()

                        open_long_seq = open(long_seqs +'/'+ long_seq, 'r')
                        read_long_seq = open_long_seq.read()

                        fix_short_seq = convert_short_seq(direction, read_manip_file)

                        trim = trim_long_seq(start, stop, read_long_seq)
                        
                        combined = trim + '\n' + fix_short_seq

                        output = open(outdir +'/'+ 'combined_'+ manip_locus + '-' + manip_taxon, 'w')
                        output.write(combined)
                        output.close()
                        open_manip_file.close()
                        open_long_seq.close() 
                        

                    



def main():
    args = parse_args()

    path_to_manip_folder = os.path.realpath(args.manipulate_seqs_folder)
    manip_folder_contents = os.listdir(path_to_manip_folder)

    path_to_long_seqs_folder = os.path.realpath(args.long_seqs_folder)
    long_seqs_folder_contents = os.listdir(path_to_long_seqs_folder)

    path_to_blast_folder = os.path.realpath(args.blast_files_folder)
    blast_folder_contents = os.listdir(path_to_blast_folder)
    assert len(blast_folder_contents) > 0

    direction_check = check_alignment_direction(blast_folder_contents, path_to_blast_folder, args.cluster_id)

    give_directions = analyze_direction(direction_check)

    find_resize = find_long_seq_size(direction_check, give_directions)
    print(find_resize)

    resize = resize_long(find_resize, give_directions, path_to_long_seqs_folder, path_to_manip_folder, args.output_dir)


if __name__ == '__main__':
    main()