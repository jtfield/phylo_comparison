#! /bin/bash
# program to reproduceably test pychorder output trees and alignments
# against trees and alignments produced by the traditional methods of inference
# Make sure all paths to Phycorder, the read file and the Useful_scripts directory are correct in the config file
# Sample command: ./comparison_commands.sh ./comparison_config.cfg


PHY_COMPARE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# get paths to the full read directory and phycorder directory
source $1

mkdir first_5_assembly

taxon_splitter.py -m --taxa_dir $master_reads --new_dir first_5_assembly --max_num 5

# run gon_phyling to assemble the first tree for phycorder to use
$phycorder_path/gon_phyling.sh $phycorder_path/comparison_5_gon_phyling.cfg

mkdir phyorder_5_plus_20

taxon_splitter.py move 25 taxa files to phycorder_5_plus_20

taxon_splitter.py delete 5 original taxa files
