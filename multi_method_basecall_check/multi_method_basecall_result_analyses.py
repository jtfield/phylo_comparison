#! /usr/bin/env python3
import argparse
import os
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
    
    # go through each methods output folder and get resulting blast results
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
    
    print(rapup_blast_results)
    read_rapup_tree = dendropy.Tree.get(path = rapup_tree, schema='newick')
    print(read_rapup_tree)











if __name__ == '__main__':
    main()
