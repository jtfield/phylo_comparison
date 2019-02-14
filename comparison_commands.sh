#! /bin/bash
# program to reproduceably test pychorder output trees and alignments
# against trees and alignments produced by the traditional methods of inference
# Make sure all paths to Phycorder, the read file and the Useful_scripts directory are correct in the config file
# Sample command: ./comparison_commands.sh ./comparison_config.cfg


PHY_COMPARE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# get paths to the full read directory and phycorder directory
source $1

# PHYCORDER RUNS BLOCK

mkdir first_5_assembly

$phycorder_path/taxon_splitter.py -m --taxa_dir $master_reads --new_dir $PHY_COMPARE/first_5_assembly --max_num 5

# # run gon_phyling to assemble the first tree for phycorder to use
$phycorder_path/gon_phyling.sh $PHY_COMPARE/comparison_5_gon_phyling.cfg
#
mkdir phycorder_5_plus_20
#
$phycorder_path/taxon_splitter.py -m --taxa_dir $master_reads --new_dir $PHY_COMPARE/phycorder_5_plus_20 --max_num 25
#
$phycorder_path/taxon_splitter.py -d --taxa_dir $PHY_COMPARE/phycorder_5_plus_20/ --max_num 5

$phycorder_path/multi_map.sh $PHY_COMPARE/phycorder_20_to_5.cfg

mkdir phycorder_25_plus_25_first

$phycorder_path/taxon_splitter.py -m --taxa_dir $master_reads --new_dir $PHY_COMPARE/phycorder_25_plus_25_first --max_num 50

$phycorder_path/taxon_splitter.py -d --taxa_dir $PHY_COMPARE/phycorder_25_plus_25_first/ --max_num 25

$phycorder_path/multi_map.sh $PHY_COMPARE/phycorder_25_plus_25_first.cfg
