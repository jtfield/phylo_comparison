#! /usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict
import numpy as np
import pandas as pd
from random import *
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import dendropy


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder')
    return parser.parse_args()

def make_df(folder_path, input_folder):
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    taxon_name = "<BlastOutput_query-def>(.+)</BlastOutput_query-def>"
    name_compile = re.compile(taxon_name)
    

    taxa_names = []
    cluster_names = []

    file_count = 0
    for file_name in input_folder:
        cluster_name_search = re.findall(cluster_compile, file_name)
        if cluster_name_search:
            if cluster_name_search[0] not in cluster_names:
                cluster_names.append(cluster_name_search[0])
        file_count+=1
        read_results = open(folder_path + "/" + file_name, "r")
        results_string = read_results.read()
        split_file = results_string.split('\n')
        name_search = re.findall(name_compile, results_string)
        if name_search[0] not in taxa_names:
            taxa_names.append(name_search[0])

    #print(taxa_names)
    #print(cluster_names)
    
    df = pd.DataFrame(columns=cluster_names, index=taxa_names)
    return df
        

        



def basecall_method_checker(folder_path, input_folder, df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    #print(df) 

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        
        # get the name of the cluster in the current file name
        # add the cluster name to the list of cluster names
        cluster_name_search = re.findall(cluster_compile, file_name)
        #tax_name_search = re.findall(name_compile, results_string)

        read_results = open(folder_path + "/" + file_name, "r")
        results_string = read_results.read()
        split_file = results_string.split('\n')

        seq_len_match = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
        taxon_name = "<BlastOutput_query-def>(.+)</BlastOutput_query-def>"
        hsp_num_1 = "<Hit_num>1</Hit_num>"
        hsp_midline = "<Hsp_midline>(.+)</Hsp_midline>"
        hit_end = "</Hsp>"
         
        hsp_compile = re.compile(hsp_num_1)
        name_compile = re.compile(taxon_name)
        len_compile = re.compile(seq_len_match)
        midline_compile = re.compile(hsp_midline)
        hit_end_compile = re.compile(hit_end)
        
        tax_name_search = re.findall(name_compile, results_string)
        #first_hit_search = re.findall(hsp_compile, line)
        midline_search = re.findall(midline_compile, results_string)
        hit_end_search = re.findall(hit_end_compile, results_string)
        
        correctly_called_bases = 0
        miscalled_bases = 0
        if cluster_name_search:
            if tax_name_search: 
                if midline_search:
                    #print(midline_search)
                    for num, midline in enumerate(midline_search[0]):
                        if midline == ' ':
                            miscalled_bases+=1
                            #miscall_base_positions["miscall_base_positons"].append(num)
                            #miscalled_base_positions.append(num)
                        elif midline == '|':
                            correctly_called_bases+=1

        df.loc[tax_name_search[0], cluster_name_search[0]] = miscalled_bases
    
    return df
   


def fig_gen(df_1, method_1, df_2, method_2, df_3, method_3):

    #hist = df.hist(bins=10)
    df_1['sums'] = df_1.sum(axis=1)
    df_2['sums'] = df_2.sum(axis=1)
    df_3['sums'] = df_3.sum(axis=1)
    sums_1 = df_1['sums']
    sums_2 = df_2['sums']
    sums_3 = df_3['sums']

    
    #plt.plot(df['sums'])
    #num_bins = 15
    #n, bins, patches = plt.hist(df['sums'], num_bins, facecolor='blue', alpha=0.5)
    #plt.show()
    #file_name = method + "_miscalled_bases.png"
    #plt.savefig(file_name)


    #stacked plots
    #fig, (ax1, ax2, ax3) = plt.subplots(1, 2)
    #fig.suptitle('Horizontally stacked subplots')
    #n, bins, patches = plt.hist(df['sums'], num_bins, facecolor='blue', alpha=0.5)
    #ax1.plot(df_1['sums'])
    #ax2.plot(df_2['sums'])
    #ax3.plot(df_3['sums'])



    fig, axs = plt.subplots(3, sharex=True, sharey=True)
    fig.suptitle('Program Miscalled Bases')
    axs[0].plot(df_1['sums'])
    axs[0].set_title('RapUp')
    
    axs[1].plot(df_2['sums'])
    axs[1].set_title('Snippy')
    
    axs[2].plot(df_3['sums'])
    axs[2].set_title('Gon_phy')
    plt.xticks([])
#    plt.show()


    #OVERLAYED DOT PLOT

#    plt.plot(sums_1, 'bo')
#    plt.plot(sums_2,'go')
#    plt.plot(sums_3, 'co')
#    plt.xticks([])
    plt.show()



    return "done"



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
   
    rapup_df = make_df(rapup_results, rapup_blast_results)
    #print(rapup_df)

    snippy_df = make_df(snippy_results, snippy_blast_results)
    #print(snippy_df)

    gon_phy_df = make_df(gon_phy_results, gon_phy_blast_results)
    #print(gon_phy_df)

    #print(rapup_blast_results)
    read_rapup_tree = dendropy.Tree.get(path = rapup_tree, schema='newick')
    
    print("rapup results") 
    rapup_basecall_check = basecall_method_checker(rapup_results, rapup_blast_results, rapup_df) 
    #print(rapup_basecall_check)
    #print(rapup_basecall_check["miscalled_bases"])
    
    print("snippy results")
    snippy_basecall_check = basecall_method_checker(snippy_results, snippy_blast_results, snippy_df)
    #print(snippy_basecall_check)
    #print(snippy_basecall_check["miscalled_bases"])

    print("gon_phy results")
    gon_phy_basecall_check = basecall_method_checker(gon_phy_results, gon_phy_blast_results, gon_phy_df)
    #print(gon_phy_basecall_check)
    #print(gon_phy_basecall_check["miscalled_bases"])

    #rapup_fig = fig_gen(rapup_basecall_check, "rapup")
    #print(rapup_fig)

    #snippy_fig = fig_gen(snippy_basecall_check, "snippy")

    #gon_phy_fig = fig_gen(gon_phy_basecall_check, "gon_phy")

    combo_fig = fig_gen(rapup_basecall_check, "rapup", snippy_basecall_check, "snippy", gon_phy_basecall_check, "gon_phy")





if __name__ == '__main__':
    main()
