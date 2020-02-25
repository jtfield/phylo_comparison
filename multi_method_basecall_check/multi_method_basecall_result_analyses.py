#! /usr/bin/env python3
import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib as plt
import seaborn as sns
import dendropy


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder')
    return parser.parse_args()


def main():
    args = parse_args()
    
    # get path to folder that contains all blast outputs for each method
    path_to_output_folder = os.path.realpath(args.output_folder)
    #print(path_to_output_folder)
    
    # go through each methods output folder and get blast result files, tree files and the reference sequence file
    rapup_results = path_to_output_folder + '/rapup_basecall/blast_results'
    snippy_results = path_to_output_folder + '/snippy_basecall/blast_results'
    gon_phy_results = path_to_output_folder + '/gon_phy_basecall/blast_results'
    
    ref_file = path_to_output_folder + '/update_alignment_dir/update_alignment_dir/alignment_ref.fas'

    rapup_blast_results = os.listdir(rapup_results)
    snippy_blast_results = os.listdir(snippy_results)
    gon_phy_blast_results = os.listdir(gon_phy_results)

    gon_phy_tree = path_to_output_folder + '/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/RAxML_bestTree.core_genome_run.out'
    rapup_tree = path_to_output_folder + '/rapup_run/combine_and_infer/RAxML_bestTree.consensusFULL'
    snippy_tree = path_to_output_folder + '' # RUN SNIPPY ALIGNMENT THROUGH RAXML
    
    #print(rapup_blast_results)
    read_rapup_tree = dendropy.Tree.get(path = rapup_tree, schema='newick')
    

    # begin adding analyzing information and sorting it based on which method and which reference was used (if applicable)
    # dendropy find the closest related taxa to the reference sequence
    #pdm = PhylogeneticDistanceMatrix.from_tree(read_rapup_tree)
    #pdm = dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix(is_store_path_edges=False)
    #pdm.compile_from_tree(read_rapup_tree)
    
    #print(pdm.as_data_table())

    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_compile = re.compile(cluster_names)
    cluster_name_list = []
    #for file_name in rapup_blast_results:
    #    cluster_name_search = re.findall(cluster_compile, file_name)
    #    if cluster_name_search:
    #        if cluster_name_search[0] not in cluster_name_list:
    #            cluster_name_list.append(cluster_name_search[0])
    #print(cluster_name_list)

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in rapup_blast_results:
        hit_and_hsp_count = 0

        cluster_name_search = re.findall(cluster_compile, file_name)
        if cluster_name_search:
            if cluster_name_search[0] not in cluster_name_list:
                cluster_name_list.append(cluster_name_search[0])
        read_results = open(rapup_results + "/" + file_name, "r")
            
        seq_len_match = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
        taxon_name = "<BlastOutput_query-def>(.+)</BlastOutput_query-def>"
        hsp_num_1 = "<Hit_num>(1)</Hit_num>"
        hsp_midline = "<Hsp_midline>(.+)</Hsp_midline>"
        hit_end = "</Hsp>"

        hsp_compile = re.compile(hsp_num_1)
        name_compile = re.compile(taxon_name)
        len_compile = re.compile(seq_len_match)
        midline_compile = re.compile(hsp_midline)
        hit_end_compile = re.compile(hit_end)
        #print(len_compile) 
        
        for line in read_results:
            name_search = re.findall(name_compile, line)
            len_search = re.findall(len_compile, line)
            first_hit_search = re.findall(hsp_compile, line)
            midline_search = re.findall(midline_compile, line)
            hit_end_search = re.findall(hit_end_compile, line)

            if name_search:
                print(name_search)
            elif len_search:
                print(len_search)
            elif first_hit_search:
                print(first_hit_search)
                hit_and_hsp_count = 1
            if midline_search and hit_and_hsp_count == 1:
                print(midline_search)
            elif hit_end_search:
                hit_and_hsp_count = 0

    print(cluster_name_list)


if __name__ == '__main__':
    main()
