#! /usr/bin/python3

import sys
import os
import glob
import dendropy
from dendropy.calculate import treecompare


# input the path to the gon_phyling folder you wish to analyze
replicate_dir = "/home/vortacs/git-repos/phylo_comparison"

os.chdir(replicate_dir)

list_of_rep_dirs = []
dir_list = os.listdir(".")
for file in dir_list:
    if "combined_output-" in file:
        list_of_rep_dirs.append(file)
    else:
        continue

for dir_name in list_of_rep_dirs:
    os.chdir(dir_name)
    rep_files = os.listdir(".")
    phycord_best_trees = []
    gon_phy_best_trees = []
    for file in rep_files:
        if "RAxML_bestTree.phycorder-" in file:
            phycord_best_trees.append(file)
        elif "RAxML_bestTree.gon_phy-" in file:
            gon_phy_best_trees.append(file)
        else:
            continue

    phycord_best_trees.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    gon_phy_best_trees.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    for phycord_count, phycord_tree in enumerate(phycord_best_trees):
        for gon_phy_count, gon_phy_tree in enumerate(gon_phy_best_trees):
            if phycord_count == gon_phy_count:

                # import the new gon_phyling produced tree
                gon_phy_tree_dp = dendropy.Tree.get(
                path = gon_phy_tree,
                schema = 'newick')

                # import the new phycorder produced tree
                phycorder_tree_dp = dendropy.Tree.get(
                path = phycord_tree,
                schema = 'newick')

                # establish common taxon namespace
                tns = dendropy.TaxonNamespace()

                print(type(gon_phy_tree))
                print(type(gon_phy_tree_dp))

                # ensure all trees loaded use common namespace
                tree1 = dendropy.Tree.get(
                data=gon_phy_tree_dp,
                schema='newick',
                taxon_namespace=tns)

                #tree2 = dendropy.Tree.get(
                #data = phycorder_tree,
                #schema = 'newick',
                #taxon_namespace = tns)

                # Unweighted Robinson-Foulds distance
                #print('\n')
                #print("UNWEIGHTED RF distance comparison between majority rule consensus rapid-updating method (phycorder)")
                #print("and traditional phylogenetics method: ")
                #print("RF: ")
                #majority_rule = treecompare.symmetric_difference(tree1, tree2)
                #print(majority_rule)

    os.chdir(replicate_dir)
