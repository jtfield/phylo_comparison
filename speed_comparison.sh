#!/bin/bash
# Phycorder and Gon_phyling timing comparison

source $1

mkdir -p $outdir

cd $outdir

# split read names for the phycorder starting tree into a file
ls "${master_reads}"/*$r1_tail | head -10 > start_tree_files.txt

# split reads for the rest of the updating using phycorder into a seperate file
ls "${master_reads}"/*$r1_tail | tail -n +11 > update_read_list.txt

start_dir=start_tree_dir
update_dir=update_alignment_dir

mkdir $start_dir
mkdir $update_dir

# link files for the starting tree into a seperate directory
for i in $(cat $outdir"/"start_tree_files.txt); do
  ln -s "${read_dir}"/$i $start_dir
  ln -s "${read_dir}"/${i%$r1_tail}$r2_tail $start_dir
done

# link files for updating the starting alignment to a seperate directory
for i in $(cat $outdir"/"update_read_list.txt); do
  ln -s "${read_dir}"/$i $update_dir
  ln -s "${read_dir}"/${i%$r1_tail}$r2_tail $update_dir
done

# create a config file on the fly with current settings
cat <<EOF > first_tree_assembly.cfg
## Config file for gon_phyling pipeline
# change the values of the variables to control the pipeline
# the variables are structured like this:
# variable_name="value_of_the_variable"
# value of the variable can be a path to a file or directory of files

# Path to the reference genome required for Parsnp
ref_genome="NONE"

# Path to the directory of reads for assembly
read_dir="$outdir/$start_dir"

# number of threads for Spades assembly and RAxML inference
threads="$THREADS"

# File stubs required to assemble paired end reads
r1_tail="$r1_tail"
r2_tail="$r2_tail"

#bootstrapping
bootstrapping="OFF"
EOF

# run gon_phyling to create the first starting tree
# time this in order to subtract the amount of time it takes for this assembly
# from the total run of gon_phyling later on
# this keeps phycorder from receiving an advantage in time immediately

# time $phycorder_path/gon_phyling.sh ./first_tree_assembly.cfg

{ time $phycorder_path/gon_phyling.sh ./first_tree_assembly.cfg 2> start_tree_make.stderr ; } 2> start_tree_make_time.txt

mkdir phycorder_required_files

cp $outdir/$start_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas $outdir/phycorder_required_files/

cp $outdir/$start_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/RAxML_bestTree.core_genome_run.out $outdir/phycorder_required_files/

if [ "$tree" != "NONE" ]; then
  cat <<phy_loop > basic.cfg

  ## Welcome to the Phycorder config file
  # Change the variable values to match the files and numbers you wish to use
  # this is the path to the phycorder directory. Dont move it.

  # PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

  # the path to the previously generated alignment file

  align="$outdir/phycorder_required_files/combo.fas"

  # the path to the previously generated tree file made from the alignment file

  tree="$outdir/phycorder_required_files/RAxML_bestTree.core_genome_run.out"

  # the directory of paired end read pairs you belonging to taxa you wish to add to your phylogeny

  read_dir="$outdir/$update_dir"

  # the number of taxa that can be added to the phylogeny at a single time
  # should be less than the number of cores available
  # should be balanced with the number of threads you will assign to the programs within phycorder

  phycorder_runs="$phycorder_runs"

  # the number of threads you wish to make available to each phycorder run
  # for mapping with bowtie and inference with RAxML

  threads="$phycorder_threads"

  # the tail identifiers of the read pairs
  # if the full read name is "Clade_1_01_R1_001.fastq" and Clade_1_01_R2_001.fastq"
  # then only put the portion of the file names that change to correspond to the read pairs
  # in this example, Clade_1_01_ identify the taxons and so must not be included

  r1_tail="$r1_tail"
  r2_tail="$r2_tail"

  # the output directory for your final information

  outdir="$outdir/phycorder-out"

  #bootstrapping
  bootstrapping=$bootstrapping

phy_loop

elif [ "$tree" == "NONE" ]; then
  cat <<phy_loop > basic.cfg

  ## Welcome to the Phycorder config file
  # Change the variable values to match the files and numbers you wish to use
  # this is the path to the phycorder directory. Dont move it.

  # PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

  # the path to the previously generated alignment file

  align="$outdir/phycorder_required_files/combo.fas"

  # the path to the previously generated tree file made from the alignment file

  tree="NONE"

  # the directory of paired end read pairs you belonging to taxa you wish to add to your phylogeny

  read_dir="$outdir/$update_dir"

  # the number of taxa that can be added to the phylogeny at a single time
  # should be less than the number of cores available
  # should be balanced with the number of threads you will assign to the programs within phycorder

  phycorder_runs="$phycorder_runs"

  # the number of threads you wish to make available to each phycorder run
  # for mapping with bowtie and inference with RAxML

  threads="$phycorder_threads"

  # the tail identifiers of the read pairs
  # if the full read name is "Clade_1_01_R1_001.fastq" and Clade_1_01_R2_001.fastq"
  # then only put the portion of the file names that change to correspond to the read pairs
  # in this example, Clade_1_01_ identify the taxons and so must not be included

  r1_tail="$r1_tail"
  r2_tail="$r2_tail"

  # the output directory for your final information

  outdir="$outdir/phycorder-out"

  #bootstrapping
  bootstrapping=$bootstrapping

phy_loop

fi

{ time $phycorder_path/multi_map.sh ./basic.cfg 2> update.stderr ; } 2> update_time.txt

# time $phycorder_path/multi_map.sh ./basic.cfg

# link files for the starting tree into the update directory for the complete amount of taxa
for i in $(cat $outdir"/"start_tree_files.txt); do
  ln -s "${read_dir}"/$i $update_dir
  ln -s "${read_dir}"/${i%$r1_tail}$r2_tail $update_dir
done

# create a config file on the fly with current settings
cat <<EOF > gon_phyling_full_assembly.cfg
## Config file for gon_phyling pipeline
# change the values of the variables to control the pipeline
# the variables are structured like this:
# variable_name="value_of_the_variable"
# value of the variable can be a path to a file or directory of files

# Path to the reference genome required for Parsnp
ref_genome="NONE"

# Path to the directory of reads for assembly
read_dir="$outdir/$update_dir"

# number of threads for Spades assembly and RAxML inference
threads="$THREADS"

# File stubs required to assemble paired end reads
r1_tail="$r1_tail"
r2_tail="$r2_tail"

#bootstrapping
bootstrapping="$bootstrapping"
EOF


{ time $phycorder_path/gon_phyling.sh ./gon_phyling_full_assembly.cfg 2> full_assembly_run.stderr ; } 2> full_assembly_time.txt

printf "Speed comparison complete. Check output files for times"
