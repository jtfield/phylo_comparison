#! /usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict
import numpy as np
import pandas as pd
from random import *
import dendropy
from dendropy.calculate import treecompare

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder')
    parser.add_argument('--prefix', default='')
    return parser.parse_args()

def make_df(folder_path, input_folder):
    which_files = folder_path.split('/')[-2]
    # print(which_files)

    gon_rap_cluster = '-(cluster\d+-cluster\d+)--'
    snip_cluster = 'combined-(cluster\d+)--'
    taxon_name = '--(.+)$'

    compile_gon_rap_cluster = re.compile(gon_rap_cluster)
    compile_snip_cluster = re.compile(snip_cluster)
    compile_taxon_name = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    result_file_count = 0
    for result_file in input_folder:
        result_file_count+=1
        # check if file has a gon_phyling to rapup cluster syntax
        find_gon_rap = re.findall(compile_gon_rap_cluster, result_file)
        # check if file has snippy to other method syntax
        find_snip = re.findall(compile_snip_cluster, result_file)
        # get the taxon name
        find_name = re.findall(compile_taxon_name, result_file)

        if find_gon_rap and find_name:
            # print(result_file)
            # print(find_name)
            # print(find_gon_rap)
            taxa_names.append(find_name[0])
            cluster_names.append(find_gon_rap[0])

        elif find_snip and find_name:
            # print(result_file)
            # print(find_name)
            # print(find_snip)
            taxa_names.append(find_name[0])
            cluster_names.append(find_snip[0])

    assert len(taxa_names) == result_file_count
    assert len(cluster_names) == result_file_count

    df = pd.DataFrame(columns=cluster_names, index=taxa_names)

    return df

def basecall_miscall_checker(folder_path, input_folder, df):
    gon_rap_cluster = '-(cluster\d+-cluster\d+)--'
    snip_cluster = 'combined-(cluster\d+)--'
    taxon_name = '--(.+)$'

    compile_gon_rap_cluster = re.compile(gon_rap_cluster)
    compile_snip_cluster = re.compile(snip_cluster)
    compile_taxon_name = re.compile(taxon_name)

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(compile_taxon_name, file_name)
        cluster_snip_search = re.findall(compile_snip_cluster, file_name)
        cluster_gon_rap_search = re.findall(compile_gon_rap_cluster, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        # results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(7)]
            #print(head)
            identical_nucs = head[2].strip('\n')
            miscalled_bases = head[4].strip('\n')
            gaps = head[6].strip('\n')
            if taxon_name_search:
                if cluster_snip_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    assert '\n' not in miscalled_bases
                    df.loc[taxon_name_search[0], cluster_snip_search[0]] = miscalled_bases
                elif cluster_gon_rap_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    assert '\n' not in miscalled_bases
                    df.loc[taxon_name_search[0], cluster_gon_rap_search[0]] = miscalled_bases

    #print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    return df


def basecall_gap_checker(folder_path, input_folder, df):
    gon_rap_cluster = '-(cluster\d+-cluster\d+)--'
    snip_cluster = 'combined-(cluster\d+)--'
    taxon_name = '--(.+)$'

    compile_gon_rap_cluster = re.compile(gon_rap_cluster)
    compile_snip_cluster = re.compile(snip_cluster)
    compile_taxon_name = re.compile(taxon_name)

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(compile_taxon_name, file_name)
        cluster_snip_search = re.findall(compile_snip_cluster, file_name)
        cluster_gon_rap_search = re.findall(compile_gon_rap_cluster, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        # results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(7)]
            #print(head)
            identical_nucs = head[2].strip('\n')
            miscalled_bases = head[4].strip('\n')
            gaps = head[6].strip('\n')
            if taxon_name_search:
                if cluster_snip_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    assert '\n' not in gaps
                    df.loc[taxon_name_search[0], cluster_snip_search[0]] = gaps
                elif cluster_gon_rap_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    assert '\n' not in gaps
                    df.loc[taxon_name_search[0], cluster_gon_rap_search[0]] = gaps

    # print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    return df


def basecall_identical_nucs_checker(folder_path, input_folder, df):
    gon_rap_cluster = '-(cluster\d+-cluster\d+)--'
    snip_cluster = 'combined-(cluster\d+)--'
    taxon_name = '--(.+)$'

    compile_gon_rap_cluster = re.compile(gon_rap_cluster)
    compile_snip_cluster = re.compile(snip_cluster)
    compile_taxon_name = re.compile(taxon_name)

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(compile_taxon_name, file_name)
        cluster_snip_search = re.findall(compile_snip_cluster, file_name)
        cluster_gon_rap_search = re.findall(compile_gon_rap_cluster, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        # results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(7)]
            #print(head)
            identical_nucs = head[2].strip('\n')
            miscalled_bases = head[4].strip('\n')
            gaps = head[6].strip('\n')
            if taxon_name_search:
                if cluster_snip_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    assert '\n' not in identical_nucs
                    df.loc[taxon_name_search[0], cluster_snip_search[0]] = identical_nucs
                elif cluster_gon_rap_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    assert '\n' not in identical_nucs
                    df.loc[taxon_name_search[0], cluster_gon_rap_search[0]] = identical_nucs

    # print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    return df

def basecall_total_nucs_checker(folder_path, input_folder, df):
    gon_rap_cluster = '-(cluster\d+-cluster\d+)--'
    snip_cluster = 'combined-(cluster\d+)--'
    taxon_name = '--(.+)$'

    compile_gon_rap_cluster = re.compile(gon_rap_cluster)
    compile_snip_cluster = re.compile(snip_cluster)
    compile_taxon_name = re.compile(taxon_name)

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(compile_taxon_name, file_name)
        cluster_snip_search = re.findall(compile_snip_cluster, file_name)
        cluster_gon_rap_search = re.findall(compile_gon_rap_cluster, file_name)

        read_results = open(folder_path + "/" + file_name, "r")
        # results_string = read_results.read()
        #split_file = results_string.split('\n')

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(9)]
            #print(head)
            identical_nucs = head[2].strip('\n')
            miscalled_bases = head[4].strip('\n')
            gaps = head[6].strip('\n')
            total_nucs = head[8].strip('\n')
            if taxon_name_search:
                if cluster_snip_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    assert '\n' not in total_nucs
                    df.loc[taxon_name_search[0], cluster_snip_search[0]] = total_nucs
                elif cluster_gon_rap_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    #add miscalled data to df under cluster and taxon name
                    assert '\n' not in total_nucs
                    df.loc[taxon_name_search[0], cluster_gon_rap_search[0]] = total_nucs

    # print(df)
    df = df.astype(int)
    df['sums'] = df.sum(axis=1)
    df['mean'] = df.mean(axis=1)
    return df


def main():
    args = parse_args()
   
    prefix = args.prefix

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
 
    gon_phy_tree = open(path_to_output_folder + '/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/fixed_gon_phy_MR.tre', 'r').read()
    rapup_tree = open(path_to_output_folder + '/rapup_run/combine_and_infer/fixed_rapup_MR.tre','r').read()
    snippy_tree = open(path_to_output_folder + '/fixed_snippy_MR.tre','r').read()
    
    rapup_to_gon_phy_miscall_df = make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)
    rapup_to_gon_phy_gap_df = make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)
    rapup_to_gon_phy_total_nucs_df = make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)
    # #print(rapup_total_nucs_df)

    snippy_to_gon_phy_miscall_df = make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    snippy_to_gon_phy_gap_df = make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    snippy_to_gon_phy_total_nucs_df = make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    #print(snippy_total_nucs_df)

    snippy_to_rapup_miscall_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)
    snippy_to_rapup_gap_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)
    snippy_to_rapup_total_nucs_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)
    #print(gon_phy_total_nucs_df)

    # rapup_to_gon_phy_basecall_check = basecall_miscall_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_miscall_df)
    # rapup_to_gon_phy_gap_check = basecall_gap_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_gap_df)
    # rapup_to_gon_phy_total_check = basecall_identical_nucs_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_total_nucs_df)

###################################################################################################
    #TREE COMPARISON
    tns = dendropy.TaxonNamespace()

    #name_grabber = '(\w+?):'
    #compile_name_grabber = re.compile(name_grabber)

    snippy_ref_name = ''
    with open(path_to_output_folder + '/core.ref.fa') as f:
        first_line = f.readline().strip()
        snippy_ref_name = first_line.strip('>')
    print(snippy_ref_name)
    ref = 'Reference'
    ref_compile = re.compile(ref)

    read_rapup_tree = dendropy.Tree.get(data = rapup_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    read_snippy_tree = dendropy.Tree.get(data = snippy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    read_gon_phy_tree = dendropy.Tree.get(data = gon_phy_tree, schema='newick', taxon_namespace=tns, preserve_underscores=True)
   
    print("rapup to gon_phyling RF results")
    rapup_to_gon_phy_phylo_compare = treecompare.symmetric_difference(read_gon_phy_tree, read_rapup_tree)
    print(rapup_to_gon_phy_phylo_compare)
    # rapup_error = calculate_error(len(rapup_names), rapup_phylo_compare)
    # print(rapup_error)

    print("snippy to gon_phyling RF results")
    snippy_to_gon_phy_phylo_compare = treecompare.symmetric_difference(read_gon_phy_tree, read_snippy_tree)
    print(snippy_to_gon_phy_phylo_compare)
    # snippy_error = calculate_error(len(snippy_names), snippy_phylo_compare)
    # print(snippy_error)

    print("snippy to rapup RF results")
    snippy_to_rapup_phylo_compare = treecompare.symmetric_difference(read_rapup_tree, read_snippy_tree)
    print(snippy_to_rapup_phylo_compare)
    # gon_phy_error = calculate_error(len(gon_phy_names), gon_phy_phylo_compare)
    # print(gon_phy_error)

####################################################################################################
    
    #BASECALL COMPARISON
    print("\n\n")
    print("rapup to gon_phy results")
    #check miscalls
    print("miscall results")
    rapup_to_gon_phy_basecall_check = basecall_miscall_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_miscall_df) 
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
    rapup_to_gon_phy_total_check = basecall_total_nucs_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_total_nucs_df)
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
    snippy_to_gon_phy_miscall_check = basecall_miscall_checker(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files, snippy_to_gon_phy_miscall_df)
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
    snippy_to_gon_phy_gap_check = basecall_gap_checker(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files, snippy_to_gon_phy_gap_df)
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
    snippy_to_gon_phy_total_nuc_check = basecall_total_nucs_checker(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files, snippy_to_gon_phy_total_nucs_df)
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
    snippy_to_rapup_basecall_check = basecall_miscall_checker(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files, snippy_to_rapup_miscall_df)
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

    snippy_to_rapup_total_nuc_check = basecall_total_nucs_checker(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files, snippy_to_rapup_total_nucs_df)
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
    

if __name__ == '__main__':
    main()