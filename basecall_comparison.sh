#!/bin/bash
# Phycorder and Gon_phyling basecall accuracy comparison

source $1

mkdir -p $outdir

cd $outdir

# split read names for the phycorder starting tree into a file
ls "${master_reads}"/*$r1_tail | head -10 > start_tree_files.txt

# split reads for the rest of the updating using phycorder into a seperate file
ls "${master_reads}"/*$r1_tail | tail -n +11 > update_read_list.txt

start_dir=start_tree_dir
update_dir=update_alignment_dir
loci_blast_indexes=loci_blast_indexes_dir

mkdir $start_dir
mkdir $update_dir
mkdir $loci_blast_indexes

# link files for the starting tree into a seperate directory
# seperate directory for both the complete run and the single addition of a taxon run
for i in $(cat $outdir"/"start_tree_files.txt); do
  ln -s "${read_dir}"/$i $start_dir
  ln -s "${read_dir}"/${i%$r1_tail}$r2_tail $start_dir
done

# link files for updating the starting alignment to a seperate directory
for i in $(cat $outdir"/"update_read_list.txt); do
  ln -s "${read_dir}"/$i $update_dir
  ln -s "${read_dir}"/${i%$r1_tail}$r2_tail $update_dir
done



####################################################################################33
# BEING BUILDING STARTING TREE FOR PHYCORDER

if [ $gon_phy_alignments == "OFF" ]; then


  # create a config file on the fly with current settings
  cat <<EOF > first_tree_assembly.cfg
  ## Config file for gon_phyling pipeline
  # change the values of the variables to control the pipeline
  # the variables are structured like this:
  # variable_name="value_of_the_variable"
  # value of the variable can be a path to a file or directory of files
  # Path to the reference genome required for Parsnp
  ref_genome="$ref_genome"
  # Path to the directory of reads for assembly
  read_dir="$outdir/$start_dir"
  # number of runs for gon_phyling to commence in parallel
  gon_phy_runs="$gon_phy_runs"
  # number of threads for Spades assembly and RAxML inference
  threads="$THREADS"
  # File stubs required to assemble paired end reads
  r1_tail="$r1_tail"
  r2_tail="$r2_tail"
  #bootstrapping
  bootstrapping="OFF"
  # type of output
  # this setting determines if you output a single, concatenated loci MSA fasta file,
  # multiple single locus MSA fasta files
  # or a single nexus file that includes the start and stop locations of each loci
  # EXAMPLE: output_type="LOCI"
  # EXAMPLE: output_type="LOCUS"
  # EXAMPLE: output_type="NEXUS" ~~~~~~NOT SUPPORTED YET~~~~~~
  output_type="$output_type"
  # output location for the positional dict file if enabled
  loci_positions="$loci_positions"
EOF

# run gon_phyling to create the first starting tree
  time $phycorder_path/gon_phyling.sh ./first_tree_assembly.cfg

  printf "First tree produced. Beginning rapid updating\n"

fi

cp $outdir/$start_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/locus_msa_files/*.fasta $outdir/$loci_blast_indexes/

