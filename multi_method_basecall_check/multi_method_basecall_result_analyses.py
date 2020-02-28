#! /usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict
import numpy as np
import pandas as pd
from random import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import dendropy


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder')
    return parser.parse_args()

def basecall_method_checker(folder_path, input_folder, output_dict):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    cluster_name = {"loci_names" : []}
    miscalls = {"miscalled_bases" : []}
    taxon_dict = {"taxon_names" : []}
    cluster_len = {"loci_len" : []}
    miscall_base_positions = {"miscall_positions" : []}

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        current_locus = ''
        tax_specific_miscalled_base_position = {}
        miscalled_base_positions = []
        hit_and_hsp_count = 0
        miscalled_bases = 0
        correctly_called_bases = 0
        cluster_name_search = re.findall(cluster_compile, file_name)
        if cluster_name_search:
            output_dict["loci_names"].append(cluster_name_search[0])
            current_locus = cluster_name_search[0]
        
        read_results = open(folder_path + "/" + file_name, "r")
        results_string = read_results.read()

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

        name_search = re.findall(name_compile, results_string)
        if name_search:
            #print(name_search)
            output_dict["taxon_names"].append(name_search[0])
        #print(taxon_dict)

        len_search = re.findall(len_compile, results_string)
        if len_search:
            #print(len_search)
            output_dict["loci_len"].append(int(len_search[0]))
        
        for line in read_results:
            #name_search = re.findall(name_compile, line)
            #len_search = re.findall(len_compile, line)
            first_hit_search = re.findall(hsp_compile, line)
            midline_search = re.findall(midline_compile, line)
            hit_end_search = re.findall(hit_end_compile, line)

            if first_hit_search:
                #print(first_hit_search)
                hit_and_hsp_count = 1
            if midline_search and hit_and_hsp_count == 1:
                #print(midline_search)
                for num, midline in enumerate(midline_search[0]):
                    if midline == ' ':
                        miscalled_bases+=1
                        #miscall_base_positions["miscall_base_positons"].append(num)
                        miscalled_base_positions.append(num)
                    elif midline == '|':
                        correctly_called_bases+=1
            #print(miscalled_bases)
            #print(correctly_called_bases)
            elif hit_end_search:
                hit_and_hsp_count = 0
        output_dict["miscalled_bases"].append(miscalled_bases)
        output_dict["miscall_positions"].append(miscalled_base_positions)
    return output_dict
    #print(cluster_name_list)
    #print(cluster_name)
    #print(miscalls)
    #print(taxon_dict)
    #print(cluster_len)
    #print(miscall_base_positions)

# THIS FUNCTION IS FOR DEV PURPOSES ONLY. I USED THIS BECAUSE THE TEST DATASET DOESNT INCLUDE ENOUGH BASES TO BE MISCALLED
# AND I NEEDED SOME MISCALLS TO TEST FIGURE MAKING
def add_fake_miscall_data_for_test(miscall_dict):
    random_num_list = [12,15,3,20,17,1,19,30,45,17,34,6]
    #print(miscall_dict["miscalled_bases"])
    
    for num, miscall_num in enumerate(miscall_dict["miscalled_bases"]):
        miscall_dict["miscalled_bases"][num] = random_num_list[num]
    #print(miscall_dict["miscalled_bases"])

    #print( miscall_dict['loci_len'])
    for num, miscall_num in enumerate(miscall_dict["miscalled_bases"]):
        #print(miscall_dict["miscall_positions"][num])
        miscall_base_count = 0
        for i in range(miscall_num):
            ran_num = randint(1, int(miscall_dict['loci_len'][num]))
            #print(ran_num)
            miscall_dict['miscall_positions'][num].append(ran_num)
            miscall_base_count+=1
    return(miscall_dict)


def fig_gen(data_dict):

    df = pd.DataFrame(data_dict, columns= ['loci_names', 'taxon_names','loci_len','miscalled_bases', 'miscall_positions'])
    #print(df)
    row_count = 0
    for index, row in df.iterrows():
        row_count+=1
        print(row)
        
        # set up the figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim(0,10)
        ax.set_ylim(0,10)

        # draw lines
        xmin = 1
        xmax = 9
        y = 5
        height = 1

        plt.hlines(y, xmin, xmax)
        plt.vlines(xmin, y - height / 2., y + height / 2.)
        plt.vlines(xmax, y - height / 2., y + height / 2.)
        
        plt.text(xmin - 0.1, y, '1', horizontalalignment='right')
        plt.text(xmax + 0.1, y, str(row['loci_len']), horizontalalignment='left')

        plt.savefig('foo.png')

        # CREATE A NEW DATAFRAME ON THE FLY TO HANDLE THE VARIABLE NUMBER OF MISCALLED BASE POSITIONS
        for num in row['miscall_positions']:
            print(num)

        #row.plot(x='loci_len',y='miscalled_bases',color='red')
    #df.plot(kind='scatter',x='loci_len',y='miscalled_bases',color='red')
    #plt.show()
    #plt.savefig('foo.png')
    #plt.savefig('foo_' + row_count '_.png')

    return(df)

def main():
    args = parse_args()
    
    cluster_name = {"loci_names" : []}
    miscalls = {"miscalled_bases" : []}
    taxon_dict = {"taxon_names" : []}
    cluster_len = {"loci_len" : []}
    miscall_base_positions = {"miscall_positions" : []}

    rapup_master_dict = {"loci_names" : [], "miscalled_bases" : [], "taxon_names" : [], "loci_len" : [], "miscall_positions" : []}
    snippy_master_dict = {"loci_names" : [], "miscalled_bases" : [], "taxon_names" : [], "loci_len" : [], "miscall_positions" : []}
    gon_phy_master_dict = {"loci_names" : [], "miscalled_bases" : [], "taxon_names" : [], "loci_len" : [], "miscall_positions" : []}

    # get path to folder that contains all blast outputs for each method
    path_to_output_folder = os.path.realpath(args.output_folder)
    
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
    
    rapup_basecall_check = basecall_method_checker(rapup_results, rapup_blast_results, rapup_master_dict) 
    #print(rapup_basecall_check)

    add_miscalls = add_fake_miscall_data_for_test(rapup_basecall_check) 
    print(add_miscalls)

    generate_figure = fig_gen(add_miscalls)
    

    # begin adding analyzing information and sorting it based on which method and which reference was used (if applicable)
    # dendropy find the closest related taxa to the reference sequence
    #pdm = PhylogeneticDistanceMatrix.from_tree(read_rapup_tree)
    #pdm = dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix(is_store_path_edges=False)
    #pdm.compile_from_tree(read_rapup_tree)
    
    #print(pdm.as_data_table())
#
#    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
#    cluster_names = 'cluster\d+'
#    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
#    cluster_compile = re.compile(cluster_names)
#    len_compile = re.compile(cluster_len)
#    cluster_name = {"loci_names" : []}
#    miscalls = {"miscalled_bases" : []}
#    taxon_dict = {"taxon_names" : []}
#    cluster_len = {"loci_len" : []}
#    miscall_base_positions = {"miscall_positions" : []}
#
#    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
#    for file_name in rapup_blast_results:
#        current_locus = ''
#        tax_specific_miscalled_base_position = {}
#        miscalled_base_positions = []
#        hit_and_hsp_count = 0
#        miscalled_bases = 0
#        correctly_called_bases = 0
#        cluster_name_search = re.findall(cluster_compile, file_name)
#        if cluster_name_search:
#            cluster_name["loci_names"].append(cluster_name_search[0])
#            current_locus = cluster_name_search[0]
#        
#        read_results = open(rapup_results + "/" + file_name, "r")
#        results_string = read_results.read()
#
#        seq_len_match = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
#        taxon_name = "<BlastOutput_query-def>(.+)</BlastOutput_query-def>"
#        hsp_num_1 = "<Hit_num>(1)</Hit_num>"
#        hsp_midline = "<Hsp_midline>(.+)</Hsp_midline>"
#        hit_end = "</Hsp>"
#
#        hsp_compile = re.compile(hsp_num_1)
#        name_compile = re.compile(taxon_name)
#        len_compile = re.compile(seq_len_match)
#        midline_compile = re.compile(hsp_midline)
#        hit_end_compile = re.compile(hit_end)
#        #print(len_compile)
#
#        name_search = re.findall(name_compile, results_string)
#        if name_search:
#            taxon_dict["taxon_names"].append(name_search[0])
#        #print(taxon_dict)
#
#        len_search = re.findall(len_compile, results_string)
#        if len_search:
#            cluster_len["loci_len"].append(len_search[0])
#        
#        for line in read_results:
#            #name_search = re.findall(name_compile, line)
#            #len_search = re.findall(len_compile, line)
#            first_hit_search = re.findall(hsp_compile, line)
#            midline_search = re.findall(midline_compile, line)
#            hit_end_search = re.findall(hit_end_compile, line)
#
#            if first_hit_search:
#                #print(first_hit_search)
#                hit_and_hsp_count = 1
#            if midline_search and hit_and_hsp_count == 1:
#                #print(midline_search)
#                for num, midline in enumerate(midline_search[0]):
#                    if midline == ' ':
#                        miscalled_bases+=1
#                        miscalled_base_positions.append(num)
#                    elif midline == '|':
#                        correctly_called_bases+=1
#            #print(miscalled_bases)
#            #print(correctly_called_bases)
#            elif hit_end_search:
#                hit_and_hsp_count = 0
#        miscall_base_positions["miscall_positions"].append(miscalled_base_positions)
#    #print(cluster_name_list)
#    print(cluster_name)
#    print(miscalls)
#    print(taxon_dict)
#    print(cluster_len)
#    print(miscall_base_positions)

if __name__ == '__main__':
    main()
