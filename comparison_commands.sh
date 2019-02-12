#! /bin/bash
# program to reproduceably test pychorder output trees and alignments
# against trees and alignments produced by the traditional methods of inference


PHY_COMPARE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

source $1
# master_reads=$1
#
# phycorder_path=$2

# useful_scripts_path=$3

mkdir first_5_assembly


$phycorder_path/gon_phyling.sh $phycorder_path/comparison_5_gon_phyling.cfg
