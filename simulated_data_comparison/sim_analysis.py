#! /usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import dendropy
from dendropy.calculate import treecompare


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder')
    parser.add_argument('--prefix', default='')
    return parser.parse_args()

def gon_rap_make_df(folder_path, input_folder):
    # cluster_names = 'cluster\d+-cluster\d+'
    cluster_names = 'cluster\d+'
    cluster_compile = re.compile(cluster_names)
    # taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    taxon_name = "--(.+)$"
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
    
    
    df = pd.DataFrame(columns=cluster_names, index=taxa_names)
    # print(df)

    return df

def snip_make_df(folder_path, input_folder):
    # cluster_names = 'cluster\d+'
    column_index = ['whole_genome']
    # cluster_compile = re.compile(cluster_names)
    # taxon_name = "basecall_results-cluster\d+-(.+)-.txt"
    taxon_name = "-cluster\d+--(.+)$"
    name_compile = re.compile(taxon_name)
    
    taxa_names = []
    cluster_names = []

    file_count = 0
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        # cluster_name_search = re.findall(cluster_compile, file_name)
        # if cluster_name_search:
            # if cluster_name_search[0] not in cluster_names:
                # cluster_names.append(cluster_name_search[0])
        if taxon_name_search:
            if taxon_name_search[0] not in taxa_names:
                taxa_names.append(taxon_name_search[0])
    
    
    df = pd.DataFrame(columns=column_index, index=taxa_names)
    # print(df)

    return df


# for checking the number of miscalls in a comparison output
def gon_rap_basecall_checker(folder_path, input_folder, miscall_df, gap_df, identical_nuc_df, total_nuc_df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    # cluster_names = 'cluster\d+-cluster\d+'
    cluster_names = 'cluster\d+'
    cluster_compile = re.compile(cluster_names)
    
    taxon_name = "--(.+)$"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        cluster_name_search = re.findall(cluster_compile, file_name)

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(9)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            total_nucs = head[8]
            # print(identical_nucs)
            # print(miscalled_bases)
            # print(gaps)
            if taxon_name_search:
                if cluster_name_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    # add miscalled data to df under cluster and taxon name
                    miscall_df.loc[taxon_name_search[0], cluster_name_search[0]] = int(miscalled_bases)
                    gap_df.loc[taxon_name_search[0], cluster_name_search[0]] = int(gaps)
                    identical_nuc_df.loc[taxon_name_search[0], cluster_name_search[0]] = int(identical_nucs)
                    total_nuc_df.loc[taxon_name_search[0], cluster_name_search[0]] = int(total_nucs)
    
    
    # df = df.astype(int)
    miscall_df['sums'] = miscall_df.sum(axis=1)
    miscall_df['mean'] = miscall_df.mean(axis=1)

    gap_df['sums'] = gap_df.sum(axis=1)
    gap_df['mean'] = gap_df.mean(axis=1)

    identical_nuc_df['sums'] = identical_nuc_df.sum(axis=1)
    identical_nuc_df['mean'] = identical_nuc_df.mean(axis=1)

    total_nuc_df['sums'] = total_nuc_df.sum(axis=1)
    total_nuc_df['mean'] = total_nuc_df.mean(axis=1)

    # print(miscall_df)
    # print(gap_df)
    # print(identical_nuc_df)
    # return df


    # for checking the number of miscalls in a comparison output
def snippy_basecall_checker(folder_path, input_folder, miscall_df, gap_df, identical_nuc_df, total_nuc_df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    # cluster_names = 'cluster\d+'
    # cluster_compile = re.compile(cluster_names)
    
    column_index = ['whole_genome']
    taxon_name = "--(.+)$"
    name_compile = re.compile(taxon_name)

    taxa_names = []
    cluster_names = []

    # BEGIN READING BLAST FILES AND RECORDING OUTPUT FOR EACH SEQUENCE
    for file_name in input_folder:
        taxon_name_search = re.findall(name_compile, file_name)
        # cluster_name_search = re.findall(cluster_compile, file_name)

        with open(folder_path + "/" + file_name) as myfile:
            head = [next(myfile) for x in range(9)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
            total_nucs = head[8]
            # print(identical_nucs)
            # print(miscalled_bases)
            # print(gaps)
            if taxon_name_search:
                # if cluster_name_search:
                    # print(file_name)
                    # print(miscalled_bases)
                    # add miscalled data to df under cluster and taxon name
                miscall_df.loc[taxon_name_search[0], column_index] = int(miscalled_bases)
                gap_df.loc[taxon_name_search[0], column_index] = int(gaps)
                identical_nuc_df.loc[taxon_name_search[0], column_index] = int(identical_nucs)
                total_nuc_df.loc[taxon_name_search[0], column_index] = int(total_nucs)
    
    
    # df = df.astype(int)
    miscall_df['sums'] = miscall_df.sum(axis=1)
    miscall_df['mean'] = miscall_df.mean(axis=1)

    gap_df['sums'] = gap_df.sum(axis=1)
    gap_df['mean'] = gap_df.mean(axis=1)

    identical_nuc_df['sums'] = identical_nuc_df.sum(axis=1)
    identical_nuc_df['mean'] = identical_nuc_df.mean(axis=1)

    total_nuc_df['sums'] = total_nuc_df.sum(axis=1)
    total_nuc_df['mean'] = total_nuc_df.mean(axis=1)

    # print(miscall_df)
    # print(gap_df)
    # print(identical_nuc_df)
    # return df

def basecall_comparison(path_to_files):


    # go through each methods output folder and get blast result files, tree files and the reference sequence file
    gon_phy_results = path_to_files + '/gon_phy_basecall/assessment_output'
    snippy_results = path_to_files + '/snippy_basecall/assessment_output'
    rapup_results = path_to_files + '/rapup_basecall/assessment_output'

    gon_phy_results_files = os.listdir(gon_phy_results)
    snippy_results_files = os.listdir(snippy_results)
    rapup_results_files = os.listdir(rapup_results)


    ##################################################################################################
    # making DFs
    gon_phy_miscall_df = gon_rap_make_df(gon_phy_results,gon_phy_results_files)
    gon_phy_gap_df = gon_rap_make_df(gon_phy_results,gon_phy_results_files)
    gon_phy_identical_df = gon_rap_make_df(gon_phy_results,gon_phy_results_files)
    gon_phy_total_nucs_df = gon_rap_make_df(gon_phy_results,gon_phy_results_files)

    snippy_miscall_df = snip_make_df(snippy_results, snippy_results_files)
    snippy_gap_df = snip_make_df(snippy_results, snippy_results_files)
    snippy_identical_df = snip_make_df(snippy_results, snippy_results_files)
    snippy_total_nucs_df = snip_make_df(snippy_results, snippy_results_files)

    rapup_miscall_df = gon_rap_make_df(rapup_results, rapup_results_files)
    rapup_gap_df = gon_rap_make_df(rapup_results, rapup_results_files)
    rapup_identical_df = gon_rap_make_df(rapup_results, rapup_results_files)
    rapup_total_nucs_df = gon_rap_make_df(rapup_results, rapup_results_files)

    ##########################################################################################################
    #BASECALL COMPARISON
    print("\n\n\n")
    print("gon_phy results")
    gon_phy_basecall_check = gon_rap_basecall_checker(gon_phy_results, gon_phy_results_files, gon_phy_miscall_df, gon_phy_gap_df, gon_phy_identical_df, gon_phy_total_nucs_df)
    # print(gon_phy_miscall_df)
    # print(gon_phy_gap_df)
    # print(gon_phy_identical_df)

    rapup_to_gon_phy_avg_miscalled = gon_phy_miscall_df['sums'].mean()
    print("average miscalls", rapup_to_gon_phy_avg_miscalled)
    rapup_to_gon_phy_miscalled_std = gon_phy_miscall_df.loc[:,"sums"].std()
    print("miscall standard deviation", rapup_to_gon_phy_miscalled_std)
    rapup_to_gon_phy_total_miscalls = gon_phy_miscall_df.loc[:,"sums"].sum()
    print("total miscalls", rapup_to_gon_phy_total_miscalls)
    #print(rapup_basecall_check['sums'])
    # gon_phy_basecall_check = gon_phy_miscall_df.rename(columns={'sums' : 'rapup_sums'})

    rapup_to_gon_phy_avg_gap = gon_phy_gap_df['sums'].mean()
    print("average gaps", rapup_to_gon_phy_avg_gap)
    rapup_to_gon_phy_gap_std = gon_phy_gap_df.loc[:,"sums"].std()
    print("gaps standard deviasion", rapup_to_gon_phy_gap_std)
    rapup_to_gon_phy_total_gaps = gon_phy_gap_df.loc[:,"sums"].sum()
    print("total gaps", rapup_to_gon_phy_total_gaps)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_gap_check = rapup_to_gon_phy_gap_check.rename(columns={'sums' : 'rapup_sums'})

    rapup_to_gon_phy_avg_identical_nuc = gon_phy_identical_df['sums'].mean()
    print("average identical nuleotides per taxon nuleotides", rapup_to_gon_phy_avg_identical_nuc)
    rapup_to_gon_phy_identical_nuc_std = gon_phy_identical_df.loc[:,"sums"].std()
    print("identical nucleotides standard deviation per taxon", rapup_to_gon_phy_identical_nuc_std)
    rapup_to_gon_phy_total_identical_nuc = gon_phy_identical_df.loc[:,"sums"].sum()
    print("identical nucleotides summed", rapup_to_gon_phy_total_identical_nuc)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_total_check = rapup_to_gon_phy_total_check.rename(columns={'sums' : 'rapup_sums'})

    #print(rapup_total_check)
    rapup_to_gon_phy_avg_total_nuc = gon_phy_total_nucs_df['sums'].mean()
    print("average all nuleotides per taxon nuleotides", rapup_to_gon_phy_avg_total_nuc)
    rapup_to_gon_phy_total_nuc_std = gon_phy_total_nucs_df.loc[:,"sums"].std()
    print("total nucleotides standard deviation per taxon", rapup_to_gon_phy_total_nuc_std)
    rapup_to_gon_phy_total_total_nuc = gon_phy_total_nucs_df.loc[:,"sums"].sum()
    print("total nucleotides summed", rapup_to_gon_phy_total_total_nuc)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_total_check = rapup_to_gon_phy_total_check.rename(columns={'sums' : 'rapup_sums'})

    #per-base results
    rapup_per_base_miscall = rapup_to_gon_phy_total_miscalls / rapup_to_gon_phy_total_total_nuc 
    rapup_per_base_gap = rapup_to_gon_phy_total_gaps / rapup_to_gon_phy_total_total_nuc
    #print(rapup_total_miscalls)
    #print(rapup_total_total_nuc)
    print("gon_phy miscalls per base")
    print(rapup_per_base_miscall)
    print("gon_phy gaps per base")
    print(rapup_per_base_gap)

    print("\n\n\n")
    print("snippy results")
    snippy_to_gon_phy_basecall_check = snippy_basecall_checker(snippy_results, snippy_results_files, snippy_miscall_df, snippy_gap_df, snippy_identical_df, snippy_total_nucs_df)
    # print(snippy_miscall_df)
    # print(snippy_gap_df)
    # print(snippy_identical_df)

    snippy_to_gon_phy_avg_miscalled = snippy_miscall_df['sums'].mean()
    print("average miscalls", snippy_to_gon_phy_avg_miscalled)
    snippy_to_gon_phy_miscalled_std = snippy_miscall_df.loc[:,"sums"].std()
    print("miscall standard deviation", snippy_to_gon_phy_miscalled_std)
    snippy_to_gon_phy_total_miscalls = snippy_miscall_df.loc[:,"sums"].sum()
    print("total miscalls", snippy_to_gon_phy_total_miscalls)
    #print(rapup_basecall_check['sums'])
    # gon_phy_basecall_check = gon_phy_miscall_df.rename(columns={'sums' : 'rapup_sums'})

    snippy_to_gon_phy_avg_gap = snippy_gap_df['sums'].mean()
    print("average gaps", snippy_to_gon_phy_avg_gap)
    snippy_to_gon_phy_gap_std = snippy_gap_df.loc[:,"sums"].std()
    print("gaps standard deviasion", snippy_to_gon_phy_gap_std)
    snippy_to_gon_phy_total_gaps = snippy_gap_df.loc[:,"sums"].sum()
    print("total gaps", snippy_to_gon_phy_total_gaps)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_gap_check = rapup_to_gon_phy_gap_check.rename(columns={'sums' : 'rapup_sums'})

    snippy_to_gon_phy_avg_identical_nuc = snippy_identical_df['sums'].mean()
    print("average identical nuleotides per taxon nuleotides", snippy_to_gon_phy_avg_identical_nuc)
    snippy_to_gon_phy_identical_nuc_std = snippy_identical_df.loc[:,"sums"].std()
    print("identical nucleotides standard deviation per taxon", snippy_to_gon_phy_identical_nuc_std)
    snippy_to_gon_phy_total_identical_nuc = snippy_identical_df.loc[:,"sums"].sum()
    print("identical nucleotides summed", snippy_to_gon_phy_total_identical_nuc)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_total_check = rapup_to_gon_phy_total_check.rename(columns={'sums' : 'rapup_sums'})

    #print(rapup_total_check)
    snippy_to_gon_phy_avg_total_nuc = snippy_total_nucs_df['sums'].mean()
    print("average all nuleotides per taxon nuleotides", snippy_to_gon_phy_avg_total_nuc)
    snippy_to_gon_phy_total_nuc_std = snippy_total_nucs_df.loc[:,"sums"].std()
    print("total nucleotides standard deviation per taxon", snippy_to_gon_phy_total_nuc_std)
    snippy_to_gon_phy_total_total_nuc = snippy_total_nucs_df.loc[:,"sums"].sum()
    print("total nucleotides summed", snippy_to_gon_phy_total_total_nuc)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_total_check = rapup_to_gon_phy_total_check.rename(columns={'sums' : 'rapup_sums'})

    #per-base results
    snippy_per_base_miscall = snippy_to_gon_phy_total_miscalls / snippy_to_gon_phy_total_total_nuc 
    snippy_per_base_gap = snippy_to_gon_phy_total_gaps / snippy_to_gon_phy_total_total_nuc
    #print(rapup_total_miscalls)
    #print(rapup_total_total_nuc)
    print("snippy miscalls per base")
    print(snippy_per_base_miscall)
    print("snippy gaps per base")
    print(snippy_per_base_gap)


    print("\n\n\n")
    print("rapup results")
    snippy_to_rapup_basecall_check = gon_rap_basecall_checker(rapup_results, rapup_results_files, rapup_miscall_df, rapup_gap_df, rapup_identical_df, rapup_total_nucs_df) 
    # print(rapup_miscall_df)
    # print(rapup_gap_df)
    # print(snippy_identical_df)

    snippy_to_rapup_avg_miscalled = rapup_miscall_df['sums'].mean()
    print("average miscalls", snippy_to_rapup_avg_miscalled)
    snippy_to_rapup_miscalled_std = rapup_miscall_df.loc[:,"sums"].std()
    print("miscall standard deviation", snippy_to_rapup_miscalled_std)
    snippy_to_rapup_total_miscalls = rapup_miscall_df.loc[:,"sums"].sum()
    print("total miscalls", snippy_to_rapup_total_miscalls)
    #print(rapup_basecall_check['sums'])
    # gon_phy_basecall_check = gon_phy_miscall_df.rename(columns={'sums' : 'rapup_sums'})

    snippy_to_rapup_avg_gap = rapup_gap_df['sums'].mean()
    print("average gaps", snippy_to_rapup_avg_gap)
    snippy_to_rapup_gap_std = rapup_gap_df.loc[:,"sums"].std()
    print("gaps standard deviasion", snippy_to_rapup_gap_std)
    snippy_to_rapup_total_gaps = rapup_gap_df.loc[:,"sums"].sum()
    print("total gaps", snippy_to_rapup_total_gaps)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_gap_check = rapup_to_gon_phy_gap_check.rename(columns={'sums' : 'rapup_sums'})

    snippy_to_rapup_avg_identical_nuc = rapup_identical_df['sums'].mean()
    print("average identical nuleotides per taxon nuleotides", snippy_to_rapup_avg_identical_nuc)
    snippy_to_rapup_identical_nuc_std = rapup_identical_df.loc[:,"sums"].std()
    print("identical nucleotides standard deviation per taxon", snippy_to_rapup_identical_nuc_std)
    snippy_to_rapup_total_identical_nuc = rapup_identical_df.loc[:,"sums"].sum()
    print("identical nucleotides summed", snippy_to_rapup_total_identical_nuc)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_total_check = rapup_to_gon_phy_total_check.rename(columns={'sums' : 'rapup_sums'})

    #print(rapup_total_check)
    snippy_to_rapup_avg_total_nuc = rapup_total_nucs_df['sums'].mean()
    print("average all nuleotides per taxon nuleotides", snippy_to_rapup_avg_total_nuc)
    snippy_to_rapup_total_nuc_std = rapup_total_nucs_df.loc[:,"sums"].std()
    print("total nucleotides standard deviation per taxon", snippy_to_rapup_total_nuc_std)
    snippy_to_rapup_total_total_nuc = rapup_total_nucs_df.loc[:,"sums"].sum()
    print("total nucleotides summed", snippy_to_rapup_total_total_nuc)
    #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_total_check = rapup_to_gon_phy_total_check.rename(columns={'sums' : 'rapup_sums'})

    #per-base results
    snip_rap_per_base_miscall = snippy_to_rapup_total_miscalls / snippy_to_rapup_total_total_nuc 
    snip_rap_per_base_gap = snippy_to_rapup_total_gaps / snippy_to_rapup_total_total_nuc
    #print(rapup_total_miscalls)
    #print(rapup_total_total_nuc)
    print("rapup miscalls per base")
    print(snip_rap_per_base_miscall)
    print("rapup gaps per base")
    print(snip_rap_per_base_gap)

    #############################################################################
    # Write full tables to csv
    gon_phy_gap_df.to_csv("gon_phy_gaps.csv", sep='\t')
    gon_phy_identical_df.to_csv("gon_phy_identical_nucs.csv", sep='\t')
    gon_phy_total_nucs_df.to_csv("gon_phy_total_nucs.csv", sep='\t')
    gon_phy_miscall_df.to_csv("gon_phy_miscalls.csv", sep='\t')

    snippy_miscall_df.to_csv("snippy_miscals.csv", sep='\t')
    snippy_gap_df.to_csv("snippy_gaps.csv", sep='\t')
    snippy_identical_df.to_csv("snippy_identical_nucs.csv", '\t')
    snippy_total_nucs_df.to_csv("snippy_total_nucs.csv", sep='\t')

    rapup_miscall_df.to_csv("rapup_miscalls.csv", sep='\t')
    rapup_gap_df.to_csv("rapup_gaps.csv", sep='\t')
    rapup_identical_df.to_csv("rapup_identical_nucs.csv", sep='\t')
    rapup_total_nucs_df.to_csv("rapup_total_nucs.csv", sep='\t')    




def phylogeny_comparison(path_to_files):
    gon_phy_tree = open(path_to_files + '/fixed_gon_phy_MR.tre', 'r').read()
    rapup_tree = open(path_to_files + '/fixed_rapup_MR.tre','r').read()
    snippy_tree = open(path_to_files + '/fixed_snippy_MR.tre','r').read()
    true_tree = open(path_to_files + '/simtree.tre','r').read()
    
    tns = dendropy.TaxonNamespace()

    read_rapup_tree = dendropy.Tree.get(data = rapup_tree, schema='newick', preserve_underscores=True, taxon_namespace=tns)
    read_snippy_tree = dendropy.Tree.get(data = snippy_tree, schema='newick', preserve_underscores=True, taxon_namespace=tns)
    read_gon_phy_tree = dendropy.Tree.get(data = gon_phy_tree, schema='newick', preserve_underscores=True, taxon_namespace=tns)
    read_true_tree = dendropy.Tree.get(data = true_tree, schema='newick', preserve_underscores=True, taxon_namespace=tns)

    print("rapup tree RF to true tree")
    print(treecompare.symmetric_difference(read_rapup_tree, read_true_tree))
    print("gon_phy tree RF to true tree")
    print(treecompare.symmetric_difference(read_gon_phy_tree, read_true_tree))
    print("snippy tree RF to true tree")
    print(treecompare.symmetric_difference(read_snippy_tree, read_true_tree))


    #tree_list = dendropy.TreeList()
    #tree_list.read(data=gon_phy_tree, schema="newick", preserve_underscores=True)
    #print(len(tree_list.taxon_namespace))
    #tree_list.read(data=rapup_tree, schema="newick", preserve_underscores=True)
    #print(len(tree_list.taxon_namespace))
    #tree_list.read(data=snippy_tree, schema="newick", preserve_underscores=True)
    #print(len(tree_list.taxon_namespace))

def fig_gen(miscalls, gaps, identical, total):
    miscall_df = pd.read_csv(miscalls, sep='\t')
    gaps_df = pd.read_csv(gaps, sep='\t')
    identical_df = pd.read_csv(identical, sep='\t')
    total_df = pd.read_csv(total, sep='\t')

    miscall_sum = miscall_df.loc[:,"sums"].sum()
    print(miscall_sum)

    gaps_sum = gaps_df.loc[:,"sums"].sum()
    print(gaps_sum)

    identical_sum = identical_df.loc[:,"sums"].sum()
    print(identical_sum)

    total_sum = total_df.loc[:,"sums"].sum()
    print(total_sum)
    
    height = [miscall_sum, gaps_sum, identical_sum]
    bars = ['miscall_sum', 'gaps_sum', 'identical_sum']
    explode = (0.1, 0, 0)  # explode 1st slice
    colors = ['gold', 'yellowgreen', 'lightcoral']

    plt.pie(height, explode=explode, labels=bars, colors=colors,
    autopct='%1.1f%%', shadow=True, startangle=140)
    
    plt.axis('equal')
    plt.show()

def main():
    args = parse_args()

    # get path to folder that contains all blast outputs for each method
    path_to_output_folder = os.path.realpath(args.output_folder)

    # basecall_comparison(path_to_output_folder)
    
    # phylogeny_comparison(path_to_output_folder)

    fig_gen("rapup_miscalls.csv", "rapup_gaps.csv", "rapup_identical_nucs.csv", "rapup_total_nucs.csv")
    fig_gen("snippy_miscalls.csv", "snippy_gaps.csv", "snippy_identical_nucs.csv", "snippy_total_nucs.csv")
    fig_gen("gon_phy_miscalls.csv", "gon_phy_gaps.csv", "gon_phy_identical_nucs.csv", "gon_phy_total_nucs.csv")





if __name__ == '__main__':
    main()
