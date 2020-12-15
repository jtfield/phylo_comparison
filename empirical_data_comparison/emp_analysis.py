#! /usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict
import numpy as np
import pandas as pd
import dendropy
from dendropy.calculate import treecompare


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder')
    parser.add_argument('--prefix', default='')
    return parser.parse_args()

def snip_rap_make_df(folder_path, input_folder):
    cluster_names = 'cluster\d+-cluster\d+'
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

def make_df(folder_path, input_folder):
    cluster_names = 'cluster\d+'
    cluster_compile = re.compile(cluster_names)
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
    
    
    df = pd.DataFrame(columns=cluster_names, index=taxa_names)
    # print(df)

    return df


# for checking the number of miscalls in a comparison output
def gon_rap_basecall_checker(folder_path, input_folder, miscall_df, gap_df, identical_nuc_df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = 'cluster\d+-cluster\d+'
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
            head = [next(myfile) for x in range(7)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
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
    
    print(miscall_df)
    print(gap_df)
    print((identical_nuc_df))
    # df = df.astype(int)
    # df['sums'] = df.sum(axis=1)
    # df['mean'] = df.mean(axis=1)
    # return df


    # for checking the number of miscalls in a comparison output
def snippy_basecall_checker(folder_path, input_folder, miscall_df, gap_df, identical_nuc_df):
    # LIST CLUSTERS/LOCI IN THIS ANALYSIS
    cluster_names = '-cluster\d+'
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
            head = [next(myfile) for x in range(7)]
            #print(head)
            identical_nucs = head[2]
            miscalled_bases = head[4]
            gaps = head[6]
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
    
    print(miscall_df)
    print(gap_df)
    print((identical_nuc_df))
    # df = df.astype(int)
    # df['sums'] = df.sum(axis=1)
    # df['mean'] = df.mean(axis=1)
    # return df



def main():
    args = parse_args()

    # get path to folder that contains all blast outputs for each method
    path_to_output_folder = os.path.realpath(args.output_folder)
    
    # go through each methods output folder and get blast result files, tree files and the reference sequence file
    rapup_to_gon_phy_results = path_to_output_folder + '/gon_phy_to_rapup/assessment_output'
    snippy_to_gon_phy_results = path_to_output_folder + '/gon_phy_to_snippy/assessment_output'
    snippy_to_rapup_results = path_to_output_folder + '/rapup_to_snippy/assessment_output'

    rapup_to_gon_phy_alignment_result_files = os.listdir(rapup_to_gon_phy_results)
    snippy_to_gon_phy_alignment_result_files = os.listdir(snippy_to_gon_phy_results)
    snippy_to_rapup_alignment_result_files = os.listdir(snippy_to_rapup_results)


    ##################################################################################################
    # making DFs
    rapup_to_gon_phy_miscall_df = snip_rap_make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)
    rapup_to_gon_phy_gap_df = snip_rap_make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)
    rapup_to_gon_phy_identical_nucs_df = snip_rap_make_df(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files)

    # snippy_to_gon_phy_miscall_df = make_df(snippy_to_gon_pntical_nucshy_results, snippy_to_gon_phy_alignment_result_files)
    # snippy_to_gon_phy_gap_df = make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)
    # snippy_to_gon_phy_total_nucs_df = make_df(snippy_to_gon_phy_results, snippy_to_gon_phy_alignment_result_files)

    # snippy_to_rapup_miscall_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)
    # snippy_to_rapup_gap_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)
    # snippy_to_rapup_total_nucs_df = make_df(snippy_to_rapup_results, snippy_to_rapup_alignment_result_files)

    ##########################################################################################################
    #BASECALL COMPARISON
    print("\n\n")
    print("rapup to gon_phy results")
    #check miscalls
    print("miscall results")
    rapup_to_gon_phy_basecall_check = gon_rap_basecall_checker(rapup_to_gon_phy_results,rapup_to_gon_phy_alignment_result_files, rapup_to_gon_phy_miscall_df, rapup_to_gon_phy_gap_df, rapup_to_gon_phy_identical_nucs_df) 
    
    
    # print(rapup_basecall_check)
    # rapup_to_gon_phy_avg_miscalled = rapup_to_gon_phy_basecall_check['sums'].mean()
    # print("average miscalls", rapup_to_gon_phy_avg_miscalled)
    # rapup_to_gon_phy_miscalled_std = rapup_to_gon_phy_basecall_check.loc[:,"sums"].std()
    # print("miscall standard deviation", rapup_to_gon_phy_miscalled_std)
    # rapup_to_gon_phy_total_miscalls = rapup_to_gon_phy_basecall_check.loc[:,"sums"].sum()
    # print("total miscalls", rapup_to_gon_phy_total_miscalls)
    # #print(rapup_basecall_check['sums'])
    # rapup_to_gon_phy_basecall_check = rapup_to_gon_phy_basecall_check.rename(columns={'sums' : 'rapup_sums'})



if __name__ == '__main__':
    main()