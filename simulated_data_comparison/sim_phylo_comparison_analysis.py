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
from dendropy.calculate import treecompare



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder')
    parser.add_argument('--prefix', default='')
    return parser.parse_args()

def get_taxa_names(tree_file):
    name_list = []
    name_grabber = '(\w+?):'
    compile_name_grabber = re.compile(name_grabber)
    
    grab_names = re.findall(compile_name_grabber, tree_file)

    if grab_names:
        for name in grab_names:
            name_list.append(name)

    return name_list


def make_df(folder_path, input_folder):
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    # taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    taxon_name = "-cluster\d+--(.+)$"
    name_compile = re.compile(taxon_name)
    
    taxa_names = []
    cluster_names = []

    file_count = 0
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        cluster_name_search = re.findall(cluster_compile, file_name)
        if cluster_name_search:
            if cluster_name_search[0] not in cluster_names:
                cluster_names.append(cluster_name_search[0])
        if taxon_name_search:
            if taxon_name_search[0] not in taxa_names:
                taxa_names.append(taxon_name_search[0])
        #file_count+=1
        #read_results = open(folder_path + "/" + file_name, "r")
        #results_string = read_results.read()
        #split_file = results_string.split('\n')
        #name_search = re.findall(name_compile, results_string)
        #if name_search[0] not in taxa_names:
        #    taxa_names.append(name_search[0])

    #print(taxa_names)
    #print(cluster_names)
    
    df = pd.DataFrame(columns=cluster_names, index=taxa_names)
    #print(df)

    return df
        
# for checking the number of identical nucleotides from a comparison output
def basecall_identical_nucs_checker(folder_path, input_folder, df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    # taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    taxon_name = "-cluster\d+--(.+)$"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        cluster_name_search = re.findall(cluster_compile, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(9)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            total_nucs = head[8]
            if taxon_name_search:
                if cluster_name_search:
                    #add miscalled data to df under cluster and taxon name
                    df.loc[taxon_name_search[0], cluster_name_search[0]] = total_nucs

    #print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    return df

# for checking the number of gaps in a comparison output
def basecall_gap_checker(folder_path, input_folder, df):
# LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    #taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    taxon_name = "-cluster\d+--(.+)$"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        cluster_name_search = re.findall(cluster_compile, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(9)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            if taxon_name_search:
                if cluster_name_search:
                    #add miscalled data to df under cluster and taxon name
                    df.loc[taxon_name_search[0], cluster_name_search[0]] = gaps
        

    #print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    return df


# for checking the number of miscalls in a comparison output
def basecall_method_checker(folder_path, input_folder, df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    #taxon_name = "basecall_results-(.+)-.txt"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        cluster_name_search = re.findall(cluster_compile, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        # print(file_name)
        results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(7)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            if taxon_name_search:
                if cluster_name_search:
                    #print(file_name)
                    #print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    df.loc[taxon_name_search[0], cluster_name_search[0]] = miscalled_bases
    
    #print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    return df
   



#####################################################################################################
# Snippy method nucleotide checks
#TODO: competely overhaul snippy analysis since it not longer recieves cluster/loci info

def snippy_make_df(folder_path, input_folder):
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    taxon_name = "basecall_results-(.+)-.txt"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    file_count = 0
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        #cluster_name_search = re.findall(cluster_compile, file_name)
        #if cluster_name_search:
        #    if cluster_name_search[0] not in cluster_names:
        #        cluster_names.append(cluster_name_search[0])
        if taxon_name_search:
            if taxon_name_search[0] not in taxa_names:
                taxa_names.append(taxon_name_search[0])
        #file_count+=1
        #read_results = open(folder_path + "/" + file_name, "r")
        #results_string = read_results.read()
        #split_file = results_string.split('\n')
        #name_search = re.findall(name_compile, results_string)
        #if name_search[0] not in taxa_names:
        #    taxa_names.append(name_search[0])

    #print(taxa_names)
    #print(cluster_names)

    df = pd.DataFrame(index=taxa_names, columns=['counts'])
    #print(df)
    return df

def snippy_basecall_method_checker(folder_path, input_folder, df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    #taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    taxon_name = "basecall_results-(.+)-.txt"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        #cluster_name_search = re.findall(cluster_compile, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(7)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            if taxon_name_search:
                #print(taxon_name_search)
                #print(file_name)
                #print(miscalled_bases)
                #add miscalled data to df under cluster and taxon name
                df.loc[taxon_name_search[0]] = int(miscalled_bases)

    #print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    #print(df)
    return df


def snippy_basecall_gap_checker(folder_path, input_folder, df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    #taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    taxon_name = "basecall_results-(.+)-.txt"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        #cluster_name_search = re.findall(cluster_compile, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(7)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            if taxon_name_search:
                #print(taxon_name_search)
                #print(file_name)
                #print(miscalled_bases)
                #add miscalled data to df under cluster and taxon name
                df.loc[taxon_name_search[0]] = int(gaps)


    #print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    #print(df)
    return df



#def snippy_basecall_gap_checker(folder_path, input_folder, df):
#    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
#    cluster_names = 'cluster\d+'
#    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
#    cluster_compile = re.compile(cluster_names)
#    len_compile = re.compile(cluster_len)
#    #taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
#    taxon_name = "basecall_results-(.+)-.txt"
#    name_compile = re.compile(taxon_name)
#
#    taxa_names = []
#    cluster_names = []
#
#    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
#    for file_name in input_folder:
#        taxon_name_search = re.findall(name_compile, file_name)
#        #cluster_name_search = re.findall(cluster_compile, file_name)
#
#        read_results = open(folder_path + "/" + file_name, "r")
#        results_string = read_results.read()
#        #split_file = results_string.split('\n')
#
#        with open(folder_path + "/" + file_name) as myfile:
#            head = [next(myfile) for x in range(7)]
#            #print(head)
#            identical_nucs = head[2]
#            miscalled_bases = head[4]
#            gaps = head[6]
#            if taxon_name_search:
#                #print(taxon_name_search)
#                #print(file_name)
#                #print(miscalled_bases)
#                #add miscalled data to df under cluster and taxon name
#                df.loc[taxon_name_search[0]] = int(gaps)
#
#    #print(df)
#    df = df.astype(int)
#    df['sums'] = df.sum(axis=1)
#    df['mean'] = df.mean(axis=1)
#    #print(df)
#    return df

def snippy_basecall_total_nucs_checker(folder_path, input_folder, df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    #taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    taxon_name = "basecall_results-(.+)-.txt"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        #cluster_name_search = re.findall(cluster_compile, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(9)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            total_nucs = head[8]
            if taxon_name_search:
                #print(taxon_name_search)
                #print(file_name)
                #print(miscalled_bases)
                #print(miscalled_bases)
                #print(gaps)
                #print("TOTAL NUCS")
                #add miscalled data to df under cluster and taxon name
                df.loc[taxon_name_search[0]] = int(total_nucs)

    #print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    #print(df)
    return df

######################################################################################################

#currently unused figure calculation section

def calculate_error(num_taxa_in_tree, rf):
    rf_max = 2 * ((num_taxa_in_tree * 2) - 3)
    error = rf / rf_max
    return error


def fig_gen(df_1, method_1, df_2, method_2, df_3, method_3):

    #hist = df.hist(bins=10)
    df_1['sums'] = df_1.sum(axis=1)
    df_2['sums'] = df_2.sum(axis=1)
    df_3['sums'] = df_3.sum(axis=1)
    sums_1 = df_1['sums']
    sums_2 = df_2['sums']
    sums_3 = df_3['sums']
    
    #print(sums_1)
    #print(sums_2)
    #print(sums_3)
    
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

#def get_snippy_ref_name(snippy_ref_fasta):

#####################################################################################################


def main():
    args = parse_args()
   
    prefix = args.prefix
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
    rapup_to_gon_phy_results = path_to_output_folder + '/gon_phy_to_rapup/assessment_output'
    snippy_to_gon_phy_results = path_to_output_folder + '/gon_phy_to_snippy/assessment_output'
    snippy_to_rapup_results = path_to_output_folder + '/rapup_to_snippy/assessment_output'
    
    ref_file = path_to_output_folder + '/update_alignment_dir/update_alignment_dir/alignment_ref.fas'

    rapup_to_gon_phy_alignment_result_files = os.listdir(rapup_to_gon_phy_results)
    snippy_to_gon_phy_alignment_result_files = os.listdir(snippy_to_gon_phy_results)
    snippy_to_rapup_alignment_result_files = os.listdir(snippy_to_rapup_results)

######################################################################3
#TREE READ (NECESSARY)
#    gon_phy_tree = open(path_to_output_folder + '/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/fixed_gon_phy_MR.tre', 'r').read()
#    rapup_tree = open(path_to_output_folder + '/rapup_run/combine_and_infer/fixed_rapup_MR.tre','r').read()
#    snippy_tree = open(path_to_output_folder + '/fixed_snippy_MR.tre','r').read()
##########################################################################

    rapup_to_gon_phy_miscall_df = make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)
    rapup_to_gon_phy_gap_df = make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)
    rapup_to_gon_phy_total_nucs_df = make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)
    #print(rapup_total_nucs_df)

    # snippy_to_gon_phy_miscall_df = snippy_make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    # snippy_to_gon_phy_gap_df = snippy_make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    # snippy_to_gon_phy_total_nucs_df = snippy_make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)

    snippy_to_gon_phy_miscall_df = make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    snippy_to_gon_phy_gap_df = make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    snippy_to_gon_phy_total_nucs_df = make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    #print(snippy_total_nucs_df)

    snippy_to_rapup_miscall_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)
    snippy_to_rapup_gap_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)
    snippy_to_rapup_total_nucs_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)
    #print(gon_phy_total_nucs_df)


###################################################################################################
#    #TREE COMPARISON
#    tns = dendropy.TaxonNamespace()
#
#    #name_grabber = '(\w+?):'
#    #compile_name_grabber = re.compile(name_grabber)
#
#    snippy_ref_name = ''
#    with open(path_to_output_folder + '/core.ref.fa') as f:
#        first_line = f.readline().strip()
#        snippy_ref_name = first_line.strip('>')
#    print(snippy_ref_name)
#    ref = 'Reference'
#    ref_compile = re.compile(ref)
#
#    read_rapup_tree = dendropy.Tree.get(data = rapup_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
#    read_snippy_tree = dendropy.Tree.get(data = snippy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
#    read_gon_phy_tree = dendropy.Tree.get(data = gon_phy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
   # read_true_tree = dendropy.Tree.get(data = true_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)

   # # str_true_tree = str(read_true_tree)
   # str_rapup_tree = str(read_rapup_tree)
   # str_snippy_tree = str(read_snippy_tree)
   # str_gon_phy_tree = str(read_gon_phy_tree)
   
   # str_snippy_tree = str_snippy_tree.replace(ref, snippy_ref_name)
  
   # # str_true_tree = str_true_tree.replace('.ref', '')
   # str_gon_phy_tree = str_gon_phy_tree.replace('.ref', '')
   # str_rapup_tree = str_rapup_tree.replace('.ref', '')
   # str_snippy_tree = str_snippy_tree.replace('.ref', '') 
   # #print(str_true_tree)

   # # true_names = get_taxa_names(str_true_tree)
   # rapup_names = get_taxa_names(str_rapup_tree)
   # snippy_names = get_taxa_names(str_snippy_tree)
   # gon_phy_names = get_taxa_names(str_gon_phy_tree)

   # #IF USING TREETOREADS AND GAVE SIMULATED READS A PREFIX OR LEFT THE PREFIX AS DEFAULT
   # # fixed_true_names = []
   # # fixed_true_tree = ''

   # str_rapup_tree = str_rapup_tree.replace(prefix, '')
   # str_snippy_tree = str_snippy_tree.replace(prefix, '')
   # str_gon_phy_tree = str_gon_phy_tree.replace(prefix, '')
   

   #print(str_true_tree)
   #print("\n")
   #print(str_rapup_tree)
   #print("\n")
   #print(str_snippy_tree)
   #print("\n")
   #print(str_gon_phy_tree)
   #print("\n")
   #assert len(true_names) == len(rapup_names) == len(snippy_names) == len(gon_phy_names)
   
   #print(true_names)
   #print(rapup_names)
   #print(snippy_names)
   #print(gon_phy_names)

   #print(len(true_names))
   #print(len(rapup_names))
   #print(len(snippy_names))
   #print(len(gon_phy_names))

   

   #prepare to make comparisons
   # fixed_true_tree = dendropy.Tree.get(data = str_true_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True, terminating_semicolon_required=False)

   # fixed_rapup_tree = dendropy.Tree.get(data = str_rapup_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True, terminating_semicolon_required=False)

   # fixed_snippy_tree = dendropy.Tree.get(data = str_snippy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True, terminating_semicolon_required=False)

   # fixed_gon_phy_tree = dendropy.Tree.get(data = str_gon_phy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True, terminating_semicolon_required=False)

   # # true_names = get_taxa_names(str_true_tree)
   # rapup_names = get_taxa_names(str_rapup_tree)
   # snippy_names = get_taxa_names(str_snippy_tree)
   # gon_phy_names = get_taxa_names(str_gon_phy_tree)

   #print(true_names)
   #print(rapup_names)
   #print(snippy_names)
   #print(gon_phy_names)

   # #PRUNE OUT NAMES NOT FOUND IN METHOD TREES (IMPORTANT IF A SUBSET OF SIMULATED DATASET WAS USED)
   # join_true_names = ''.join(true_names)
   # #print(join_true_names)
   # names_not_shared_list = []
   # #for name in true_names:
   # for name in rapup_names:
   #     compile_name = re.compile(name)
   #     find_name = re.findall(name, join_true_names)
   #     #print(name)
   #     assert name in snippy_names
   #     assert name in gon_phy_names
   #     if not find_name:
   #         names_not_shared_list.append(name)

   # if len(names_not_shared_list) > 0:
   #     #read_true_tree.prune_taxa_with_labels(names_not_shared_list)
   #     fixed_true_tree.prune_taxa_with_labels(names_not_shared_list)

   # print(names_not_shared_list)

   # assert len(fixed_true_tree.leaf_nodes()) == len(fixed_rapup_tree.leaf_nodes())
   # read_true_tree.write(path= path_to_output_folder + "/true_tree_subset.tre", schema="newick")
##########################################################################################################
#PIMARY COMPARE SECTION
#   print("rapup to gon_phyling RF results")
#   rapup_to_gon_phy_phylo_compare = treecompare.symmetric_difference(read_gon_phy_tree, read_rapup_tree)
#   print(rapup_to_gon_phy_phylo_compare)
#   # rapup_error = calculate_error(len(rapup_names), rapup_phylo_compare)
#   # print(rapup_error)
#
#   print("snippy to gon_phyling RF results")
#   snippy_to_gon_phy_phylo_compare = treecompare.symmetric_difference(read_gon_phy_tree, read_snippy_tree)
#   print(snippy_to_gon_phy_phylo_compare)
#   # snippy_error = calculate_error(len(snippy_names), snippy_phylo_compare)
#   # print(snippy_error)
#
#   print("snippy to rapup RF results")
#   snippy_to_rapup_phylo_compare = treecompare.symmetric_difference(read_rapup_tree, read_snippy_tree)
#   print(snippy_to_rapup_phylo_compare)
#   # gon_phy_error = calculate_error(len(gon_phy_names), gon_phy_phylo_compare)
#   # print(gon_phy_error)

###################################################################################################
    
    #BASECALL COMPARISON
    print("\n\n")
    print("rapup to gon_phy results")
    #check miscalls
    print("miscall results")
    rapup_to_gon_phy_basecall_check = basecall_method_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_miscall_df) 
    #print(rapup_basecall_check)
    rapup_to_gon_phy_avg_miscalled = rapup_to_gon_phy_basecall_check['sums'].mean()
    print("average miscalls", rapup_to_gon_phy_avg_miscalled)
    rapup_to_gon_phy_miscalled_std = rapup_to_gon_phy_basecall_check.loc[:,"sums"].std()
    print("miscall standard deviation", rapup_to_gon_phy_miscalled_std)
    rapup_to_gon_phy_total_miscalls = rapup_to_gon_phy_basecall_check.loc[:,"sums"].sum()
    print("total miscalls", rapup_to_gon_phy_total_miscalls)
    #print(rapup_basecall_check['sums'])
    rapup_to_gon_phy_basecall_check = rapup_to_gon_phy_basecall_check.rename(columns={'sums' : 'rapup_sums'})
 
    #check gaps
    rapup_to_gon_phy_gap_check = basecall_gap_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_gap_df)
    #print(rapup_gap_check)
    rapup_to_gon_phy_avg_gap = rapup_to_gon_phy_gap_check['sums'].mean()
    print("average gaps", rapup_to_gon_phy_avg_gap)
    rapup_to_gon_phy_gap_std = rapup_to_gon_phy_gap_check.loc[:,"sums"].std()
    print("gaps standard deviasion", rapup_to_gon_phy_gap_std)
    rapup_to_gon_phy_total_gaps = rapup_to_gon_phy_gap_check.loc[:,"sums"].sum()
    print("total gaps", rapup_to_gon_phy_total_gaps)
    #print(rapup_basecall_check['sums'])
    rapup_to_gon_phy_gap_check = rapup_to_gon_phy_gap_check.rename(columns={'sums' : 'rapup_sums'})
 
    #check total nucleotides
    rapup_to_gon_phy_total_check = basecall_identical_nucs_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_total_nucs_df)
    #print(rapup_total_check)
    rapup_to_gon_phy_avg_total_nuc = rapup_to_gon_phy_total_check['sums'].mean()
    print("average all nuleotides per taxon nuleotides", rapup_to_gon_phy_avg_total_nuc)
    rapup_to_gon_phy_total_nuc_std = rapup_to_gon_phy_total_check.loc[:,"sums"].std()
    print("total nucleotides standard deviation per taxon", rapup_to_gon_phy_total_nuc_std)
    rapup_to_gon_phy_total_total_nuc = rapup_to_gon_phy_total_check.loc[:,"sums"].sum()
    print("total nucleotides summed", rapup_to_gon_phy_total_total_nuc)
    #print(rapup_basecall_check['sums'])
    rapup_to_gon_phy_total_check = rapup_to_gon_phy_total_check.rename(columns={'sums' : 'rapup_sums'})
     
    #per-base results
    rapup_per_base_miscall = rapup_to_gon_phy_total_miscalls / rapup_to_gon_phy_total_total_nuc 
    rapup_per_base_gap = rapup_to_gon_phy_total_gaps / rapup_to_gon_phy_total_total_nuc
    #print(rapup_total_miscalls)
    #print(rapup_total_total_nuc)
    print("rapup miscalls per base")
    print(rapup_per_base_miscall)
    print("rapup gaps per base")
    print(rapup_per_base_gap)


    print("\n\n")
    print("snippy to gon_phyling results")
    #miscall check
    snippy_to_gon_phy_miscall_check = snippy_basecall_method_checker(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files, snippy_to_gon_phy_miscall_df)
    #print(snippy_basecall_check)
    snippy_to_gon_phy_avg_miscalled = snippy_to_gon_phy_miscall_check['sums'].mean()
    print("average miscalls per taxon", snippy_to_gon_phy_avg_miscalled)
    snippy_to_gon_phy_miscall_std = snippy_to_gon_phy_miscall_check.loc[:,"sums"].std()
    print("standard deviation of miscalls per taxon", snippy_to_gon_phy_miscall_std)
    snippy_to_gon_phy_total_miscalls = snippy_to_gon_phy_miscall_check.loc[:,"sums"].sum()
    print("summed total miscalls", snippy_to_gon_phy_total_miscalls)
    #print(snippy_basecall_check['sums'])
    snippy_basecall_check = snippy_to_gon_phy_miscall_check.rename(columns={'sums' : 'snippy_sums'})

    #gap check
    snippy_to_gon_phy_gap_check = snippy_basecall_gap_checker(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files, snippy_to_gon_phy_gap_df)
    #print(snippy_gap_check)
    snippy_to_gon_phy_avg_gap = snippy_to_gon_phy_gap_check['sums'].mean()
    print("average gaps per taxon", snippy_to_gon_phy_avg_gap)
    snippy_to_gon_phy_gap_std = snippy_to_gon_phy_gap_check.loc[:,"sums"].std()
    print("standard deviation of gaps per taxon", snippy_to_gon_phy_gap_std)
    snippy_to_gon_phy_total_gaps = snippy_to_gon_phy_gap_check.loc[:,"sums"].sum()
    print("summed total gaps", snippy_to_gon_phy_total_gaps)
    #print(snippy_basecall_check['sums'])
    snippy_to_gon_phy_gap_check = snippy_to_gon_phy_gap_check.rename(columns={'sums' : 'snippy_sums'})

    #total nucleotide check
    snippy_to_gon_phy_total_nuc_check = snippy_basecall_total_nucs_checker(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files, snippy_to_gon_phy_total_nucs_df)
    #print(snippy_basecall_check)
    snippy_to_gon_phy_avg_total_nuc = snippy_to_gon_phy_total_nuc_check['sums'].mean()
    print("average total nucleotides per taxon", snippy_to_gon_phy_avg_total_nuc)
    snippy_to_gon_phy_total_nuc_std = snippy_to_gon_phy_total_nuc_check.loc[:,"sums"].std()
    print("standard deviation of total nucleotides per taxon", snippy_to_gon_phy_total_nuc_std)
    snippy_to_gon_phy_total_total_nuc = snippy_to_gon_phy_total_nuc_check.loc[:,"sums"].sum()
    print("summed total nucleotides", snippy_to_gon_phy_total_total_nuc)
    #print(snippy_basecall_check['sums'])
    snippy_to_gon_phy_total_nuc_check = snippy_to_gon_phy_total_nuc_check.rename(columns={'sums' : 'snippy_sums'})

    #get per-base miscall and gap rate
    snippy_to_gon_phy_per_base_miscall = snippy_to_gon_phy_total_miscalls / snippy_to_gon_phy_total_total_nuc
    snippy_to_gon_phy_per_base_gap = snippy_to_gon_phy_total_gaps / snippy_to_gon_phy_total_total_nuc
    print("snippy miscalls per base")
    print(snippy_to_gon_phy_per_base_miscall)
    print("snippy gaps per base")
    print(snippy_to_gon_phy_per_base_gap)


    print("\n\n")
    print("snippy to rapup comparison results")
    snippy_to_rapup_basecall_check = basecall_method_checker(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files, snippy_to_rapup_miscall_df)
    #print(gon_phy_basecall_check)
    snippy_to_rapup_avg_miscalled = snippy_to_rapup_basecall_check['sums'].mean()
    print("average miscalls per taxon", snippy_to_rapup_avg_miscalled)
    snippy_to_rapup_miscall_std = snippy_to_rapup_basecall_check.loc[:,"sums"].std()
    print("standard deviation of miscalls per taxon", snippy_to_rapup_miscall_std)
    snippy_to_rapup_total_miscalls = snippy_to_rapup_basecall_check.loc[:,"sums"].sum()
    print("total miscalls", snippy_to_rapup_total_miscalls)
    #print(gon_phy_basecall_check['sums'])
    snippy_to_rapup_basecall_check = snippy_to_rapup_basecall_check.rename(columns={'sums' : 'gon_phy_sums'})

    snippy_to_rapup_gap_check = basecall_gap_checker(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files, snippy_to_rapup_gap_df)
    #print(gon_phy_basecall_check)
    snippy_to_rapup_avg_gap = snippy_to_rapup_gap_check['sums'].mean()
    print("average gaps per taxon", snippy_to_rapup_avg_gap)
    snippy_to_rapup_gap_std = snippy_to_rapup_gap_check.loc[:,"sums"].std()
    print("standard deviation of gaps per taxon", snippy_to_rapup_gap_std)
    snippy_to_rapup_total_gaps = snippy_to_rapup_gap_check.loc[:,"sums"].sum()
    print("summed total gaps", snippy_to_rapup_total_gaps)
    #print(gon_phy_basecall_check['sums'])
    snippy_to_rapup_gap_check = snippy_to_rapup_gap_check.rename(columns={'sums' : 'gon_phy_sums'})

    snippy_to_rapup_total_nuc_check = basecall_identical_nucs_checker(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files, snippy_to_rapup_total_nucs_df)
    #print(gon_phy_basecall_check)
    snippy_to_rapup_avg_total_nuc = snippy_to_rapup_total_nuc_check['sums'].mean()
    print("average total nucleotides per taxon", snippy_to_rapup_avg_total_nuc)
    snippy_to_rapup_total_nuc_std = snippy_to_rapup_total_nuc_check.loc[:,"sums"].std()
    print("standard deviation of total nucleotides per taxon", snippy_to_rapup_total_nuc_std)
    snippy_to_rapup_total_total_nuc = snippy_to_rapup_total_nuc_check.loc[:,"sums"].sum()
    print("summed total of nucleotides", snippy_to_rapup_total_total_nuc)
    #print(gon_phy_basecall_check['sums'])
    snippy_to_rapup_total_nuc_check = snippy_to_rapup_total_nuc_check.rename(columns={'sums' : 'gon_phy_sums'})
   
    #per-base analysis results
    snippy_to_rapup_per_base_miscall = snippy_to_rapup_total_miscalls / snippy_to_rapup_total_total_nuc
    snippy_to_rapup_per_base_gap = snippy_to_rapup_total_gaps / snippy_to_rapup_total_total_nuc
    print("gon_phy miscalls per base")
    print(snippy_to_rapup_per_base_miscall)
    print("gon_phy gaps per base")
    print(snippy_to_rapup_per_base_gap)
    
    #print(gon_phy_basecall_check['gon_phy_sums'])    
    #combine_sums_df = pd.concat([rapup_basecall_check['rapup_sums'], snippy_basecall_check['snippy_sums'], gon_phy_basecall_check['gon_phy_sums']] , axis=1, ignore_index=False, sort=True)

    #combine_index = combine_sums_df['taxa'] = combine_sums_df.index
    
    

    #rapup_fig = fig_gen(rapup_basecall_check, "rapup")
    #print(rapup_fig)

    #snippy_fig = fig_gen(snippy_basecall_check, "snippy")

    #gon_phy_fig = fig_gen(gon_phy_basecall_check, "gon_phy")

    #combo_fig = fig_gen(rapup_basecall_check, "rapup", snippy_basecall_check, "snippy", gon_phy_basecall_check, "gon_phy")





if __name__ == '__main__':
    main()
