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
    taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
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
        

        



def basecall_method_checker(folder_path, input_folder, df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+'
    cluster_len = "<BlastOutput_query-len>(\d+)</BlastOutput_query-len>"
    cluster_compile = re.compile(cluster_names)
    len_compile = re.compile(cluster_len)
    taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
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
            head = [next(myfile) for x in range(3)]
            #print(head)
            miscalled_bases = head[2]
            if taxon_name_search:
                if cluster_name_search:
                    #add miscalled data to df under cluster and taxon name
                    df.loc[taxon_name_search[0], cluster_name_search[0]] = miscalled_bases
    
    #print(df) 
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

    gon_phy_tree = open(path_to_output_folder + '/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/RAxML_bestTree.core_genome_run.out', 'r').read()
    rapup_tree = open(path_to_output_folder + '/rapup_run/combine_and_infer/RAxML_bestTree.consensusFULL','r').read()
    snippy_tree = open(path_to_output_folder + '/RAxML_bestTree.snippy_tree','r').read()
    true_tree = open(path_to_output_folder + '/true_tree.tre','r').read()

    rapup_df = make_df(rapup_results, rapup_blast_results)
    #print(rapup_df)

    snippy_df = make_df(snippy_results, snippy_blast_results)
    #print(snippy_df)

    gon_phy_df = make_df(gon_phy_results, gon_phy_blast_results)
    #print(gon_phy_df)

    #TREE COMPARISON
    #tns = dendropy.TaxonNamespace()

    #name_grabber = '(\w+?):'
    #compile_name_grabber = re.compile(name_grabber)

    #snippy_ref_name = ''
    #with open(path_to_output_folder + '/core.ref.fa') as f:
    #    first_line = f.readline().strip()
    #    snippy_ref_name = first_line.strip('>')
    #print(snippy_ref_name)

    #ref = 'Reference'
    #ref_compile = re.compile(ref)
    

    #read_rapup_tree = dendropy.Tree.get(data = rapup_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    #read_snippy_tree = dendropy.Tree.get(data = snippy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    #read_gon_phy_tree = dendropy.Tree.get(data = gon_phy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    #read_true_tree = dendropy.Tree.get(data = true_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    
    #str_true_tree = str(read_true_tree)
    #str_rapup_tree = str(read_rapup_tree)
    #str_snippy_tree = str(read_snippy_tree)
    #str_gon_phy_tree = str(read_gon_phy_tree)

    #str_snippy_tree = str_snippy_tree.replace(ref, snippy_ref_name)
    #str_gon_phy_tree = str_gon_phy_tree.replace('.ref', '')
    #str_rapup_tree = str_rapup_tree.replace('.ref', '')
    #str_snippy_tree = str_snippy_tree.replace('.ref', '')

    #true_names = get_taxa_names(str_true_tree)
    #rapup_names = get_taxa_names(str_rapup_tree)
    #snippy_names = get_taxa_names(str_snippy_tree)
    #gon_phy_names = get_taxa_names(str_gon_phy_tree)
   
    #print(true_names)
    #print(rapup_names)

    #print(str_rapup_tree)
    #print(str_gon_phy_tree)
    #print(str_snippy_tree)

    #print(len(rapup_names))
    #print(len(snippy_names))
    #print(len(gon_phy_names))

    #join_true_names = ''.join(true_names)
    #print(join_true_names)
    #names_not_shared_list = []
    #for name in rapup_names:
    #for name in true_names:
        #compile_name = re.compile(name)
        #find_name = re.findall(name, join_true_names)
        #print(name)
        #assert name in snippy_names
        #assert name in gon_phy_names
        #if name not in true_names:
    #    if name not in rapup_names:
        #if not find_name:
    #        names_not_shared_list.append(name)
    
    #print(names_not_in_true_tree)
    #if len(names_not_shared_list) > 0:
    #    read_true_tree.prune_taxa_with_labels(names_not_shared_list)

    #print(str_snippy_tree)
    #fixed_snippy_tree = dendropy.Tree.get(data = str_snippy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True, terminating_semicolon_required=False)

    #fixed_gon_phy_tree = dendropy.Tree.get(data = str_gon_phy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True, terminating_semicolon_required=False)

    #print(read_true_tree.leaf_nodes())
    #print(read_rapup_tree.leaf_nodes())
    #print(names_not_in_true_tree)

    #assert len(read_true_tree.leaf_nodes()) == len(read_rapup_tree.leaf_nodes())
    #read_true_tree.write(path= path_to_output_folder + "/true_tree_subset.tre", schema="newick")


    #BASECALL COMPARISON
    print("rapup results")
    #rapup_phylo_compare = treecompare.symmetric_difference(read_true_tree, read_rapup_tree)
    #print(rapup_phylo_compare)
    rapup_basecall_check = basecall_method_checker(rapup_results, rapup_blast_results, rapup_df) 
    #print(rapup_basecall_check)
    #print(rapup_basecall_check["miscalled_bases"])
    
    print("snippy results")
    #snippy_phylo_compare = treecompare.symmetric_difference(read_true_tree, fixed_snippy_tree)
    #print(snippy_phylo_compare)
    snippy_basecall_check = basecall_method_checker(snippy_results, snippy_blast_results, snippy_df)
    #print(snippy_basecall_check)
    #print(snippy_basecall_check["miscalled_bases"])

    print("gon_phy results")
    #gon_phy_phylo_compare = treecompare.symmetric_difference(read_true_tree, fixed_gon_phy_tree)
    #print(gon_phy_phylo_compare)
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
