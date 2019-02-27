#! /bin/bash
# program to reproduceably test pychorder output trees and alignments
# against trees and alignments produced by the traditional methods of inference
# Make sure all paths to Phycorder, the read file and the Useful_scripts directory are correct in the config file
# Sample command: ./comparison_commands.sh ./comparison_config.cfg


PHY_COMPARE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

source $1

REP_TEMPFILE=/tmp/rep_temp.tmp
echo 0 > $REP_TEMPFILE

for ((n=0;n<$REPLICATES;n++)); do

  REP_COUNTER=$[$(cat $REP_TEMPFILE) + 1]

  mkdir -p $outdir

  maind=$(pwd)

  cd $outdir

  workd=$(pwd)

  ls ${master_reads}/*$r1_tail | split -a 20 -l $file_split_num

  # potential randomization command
  #ls ${master_reads}/*$r1_tail | sort -R | split -a 20 -l $file_split_num

  wait

  for j in $(ls xa*); do
    mkdir $j-read-dir
    for i in $(cat $j); do
      cp ${read_dir}/$i $j-read-dir/
      cp ${read_dir}/${i%$r1_tail}$r2_tail $j-read-dir/
      wait
    done
    wait
  done

  first_tree_dir=$(ls -d */ | head -1)

  second_tree_dir=$(ls -d */ | head -2 | tail -1)

  printf "$first_tree_dir"

  printf "$second_tree_dir"

  wait

  cat <<EOF > first_tree_assembly.cfg
  ## Config file for gon_phyling pipeline
  # change the values of the variables to control the pipeline
  # the variables are structured like this:
  # variable_name="value_of_the_variable"
  # value of the variable can be a path to a file or directory of files

  # Path to the reference genome required for Parsnp
  ref_genome="/path/to/reference/genome"

  # Path to the directory of reads for assembly
  read_dir="$workd/$first_tree_dir"

  # number of threads for Spades assembly and RAxML inference
  threads="4"

  # File stubs required to assemble paired end reads
  r1_tail="R1.fastq"
  r2_tail="R2.fastq"
EOF

  $phycorder_path/gon_phyling.sh ./first_tree_assembly.cfg

  mkdir updated_phycorder_required_files

  cp $workd/$first_tree_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas ./updated_phycorder_required_files/

  cp $workd/$first_tree_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/RAxML_bestTree.core_genome_run.out ./updated_phycorder_required_files/

  wait

  for dir in $(ls -d */ ); do
    # printf $dir
    # printf $first_tree_dir
    if [ "$dir" == "$first_tree_dir" ]; then

      printf "skipping first first tree file"

    elif [ "$dir" == "$second_tree_dir" ]; then

      # cat <<phy_loop > $dir_phycorder_run.cfg
      cat <<phy_loop > basic.cfg

      ## Welcome to the Phycorder config file
      # Change the variable values to match the files and numbers you wish to use
      # this is the path to the phycorder directory. Dont move it.

      # PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

      # the path to the previously generated alignment file

      align="$workd/updated_phycorder_required_files/combo.fas"

      # the path to the previously generated tree file made from the alignment file

      tree="$workd/updated_phycorder_required_files/RAxML_bestTree.core_genome_run.out"

      # the directory of paired end read pairs you belonging to taxa you wish to add to your phylogeny

      read_dir="$workd/$dir"

      # the number of taxa that can be added to the phylogeny at a single time
      # should be less than the number of cores available
      # should be balanced with the number of threads you will assign to the programs within phycorder

      phycorder_runs="2"

      # the number of threads you wish to make available to each phycorder run
      # for mapping with bowtie and inference with RAxML

      threads="2"

      # the tail identifiers of the read pairs
      # if the full read name is "Clade_1_01_R1_001.fastq" and Clade_1_01_R2_001.fastq"
      # then only put the portion of the file names that change to correspond to the read pairs
      # in this example, Clade_1_01_ identify the taxons and so must not be included

      r1_tail="R1.fastq"
      r2_tail="R2.fastq"

      # the output directory for your final information

      outdir="$workd/$dir/phycorder-out"

phy_loop

      wait

      time $phycorder_path/multi_map.sh ./basic.cfg

      cp $workd/$dir/phycorder-out/combine_and_infer/RAxML_bestTree.consensusFULL ./updated_phycorder_required_files/

      cp $workd/$dir/phycorder-out/combine_and_infer/extended.aln ./updated_phycorder_required_files/

      wait
    else

      cat <<phy_loop > basic.cfg

      ## Welcome to the Phycorder config file
      # Change the variable values to match the files and numbers you wish to use
      # this is the path to the phycorder directory. Dont move it.

      # PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

      # the path to the previously generated alignment file

      align="$workd/updated_phycorder_required_files/extended.aln"

      # the path to the previously generated tree file made from the alignment file

      tree="$workd/updated_phycorder_required_files/RAxML_bestTree.consensusFULL"

      # the directory of paired end read pairs you belonging to taxa you wish to add to your phylogeny

      read_dir="$workd/$dir"

      # the number of taxa that can be added to the phylogeny at a single time
      # should be less than the number of cores available
      # should be balanced with the number of threads you will assign to the programs within phycorder

      phycorder_runs="2"

      # the number of threads you wish to make available to each phycorder run
      # for mapping with bowtie and inference with RAxML

      threads="2"

      # the tail identifiers of the read pairs
      # if the full read name is "Clade_1_01_R1_001.fastq" and Clade_1_01_R2_001.fastq"
      # then only put the portion of the file names that change to correspond to the read pairs
      # in this example, Clade_1_01_ identify the taxons and so must not be included

      r1_tail="R1.fastq"
      r2_tail="R2.fastq"

      # the output directory for your final information

      outdir="$workd/$dir/phycorder-out"

phy_loop

      wait

      time $phycorder_path/multi_map.sh ./basic.cfg

      rm ./updated_phycorder_required_files/RAxML_bestTree.consensusFULL

      rm ./updated_phycorder_required_files/extended.aln

      cp $workd/$dir/phycorder-out/combine_and_infer/RAxML_bestTree.consensusFULL ./updated_phycorder_required_files/

      cp $workd/$dir/phycorder-out/combine_and_infer/extended.aln ./updated_phycorder_required_files/

      wait
    fi

  done
  # END OF PHYCORDER RUN


  rm -r ./updated_phycorder_required_files/

  mkdir $maind/phycorder_results


  TEMPFILE=/tmp/$$.tmp
  echo 0 > $TEMPFILE

  for dir in $(ls -d */ ); do
    COUNTER=$[$(cat $TEMPFILE) + 1]

    if [ "$dir" == "$first_tree_dir" ]; then

      printf "skipping first first tree file"

      echo "$COUNTER"

      echo $COUNTER > $TEMPFILE

    else


      cp $workd/$dir/phycorder-out/combine_and_infer/RAxML_bestTree.consensusFULL $maind/phycorder_results/RAxML_bestTree.phycorder-$COUNTER-.out

      cp $workd/$dir/phycorder-out/combine_and_infer/extended.aln $maind/phycorder_results/phycorder-$COUNTER-.aln

      cp $workd/$dir/phycorder-out/combine_and_infer/RAxML_bipartitions.majority_rule_bootstrap_consensus $maind/phycorder_results/RAxML_bipartitions.phycorder_majority_rule-$COUNTER-.out

      echo "$COUNTER"

      echo $COUNTER > $TEMPFILE


    fi
    wait

  done

  unlink $TEMPFILE

  # GON PHYLING SECTION


  cd $maind

  gon_results=gon_phy_results

  gon_runs=gon_phy_runs_dir

  mkdir $gon_results

  mkdir $gon_runs

  TEMPFILE=/tmp/$$.tmp
  echo 0 > $TEMPFILE

  cd $workd
  for dir in $(ls -d */ ); do
    # if [ "$dir" != "$gon_results" ] || [ "$dir" != "$gon_runs" ]; then
    if [ "$dir" == "$first_tree_dir" ]; then

      printf "skipping first tree file"

        # Fetch the value and increase it
      COUNTER=$[$(cat $TEMPFILE) + 1]

      cp $dir/*.fastq $maind/gon_phy_runs_dir/

      echo $COUNTER > $TEMPFILE
  #     fi
  #
  #   elif [ "$dir" != "$gon_results" ] || [ "$dir" != "$gon_runs" ]; then

  else
    printf "$dir"
    # Fetch the value and increase it
    COUNTER=$[$(cat $TEMPFILE) + 1]

    cp $dir/*.fastq $maind/gon_phy_runs_dir/

    printf "copied $dir files to gon_phy_runs_dir"

    wait

    cat <<gon_phy_loop > gon_phy_basic.cfg

    ## Config file for gon_phyling pipeline
    # change the values of the variables to control the pipeline
    # the variables are structured like this:
    # variable_name="value_of_the_variable"
    # value of the variable can be a path to a file or directory of files

    # Path to the reference genome required for Parsnp
    ref_genome="/path/to/reference/genome"

    # Path to the directory of reads for assembly
    read_dir="$maind/gon_phy_runs_dir"

    # number of threads for Spades assembly and RAxML inference
    threads="4"

    # File stubs required to assemble paired end reads
    r1_tail="R1.fastq"
    r2_tail="R2.fastq"

gon_phy_loop

  #   else
  #     printf "$dir"
   $phycorder_path/gon_phyling.sh ./gon_phy_basic.cfg

    wait

    mv $maind/gon_phy_runs_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/RAxML_bestTree.core_genome_run.out $maind/gon_phy_results/RAxML_bestTree.gon_phy-$COUNTER-.out

    mv $maind/gon_phy_runs_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas $maind/gon_phy_results/gon_phy-$COUNTER-.fas

    mv $maind/gon_phy_runs_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/RAxML_bipartitions.core_genome_run.out $maind/gon_phy_results/RAxML_bipartitions.gon_phy_majority_rule-$COUNTER-.out

    echo "$COUNTER"

    echo $COUNTER > $TEMPFILE

    rm -r $maind/gon_phy_runs_dir/trimmed*

    wait

    fi

  done

  unlink $TEMPFILE

  mkdir $maind/combined_outputs

  cp $maind/gon_phy_results/* $maind/combined_outputs/

  cp $maind/phycorder_results/* $maind/combined_outputs/

  wait

  # ANALYSES SUMMARY SECTION

  # produce data outputs for gon_phyling alignments

  cd $maind/combined_outputs

  TEMPFILE=/tmp/$$.tmp
  echo 1 > $TEMPFILE
  for seq_file in $(ls *.fas); do
    COUNTER=$[$(cat $TEMPFILE) + 1]

    $phycorder_path/phy_stats.py --align_file $seq_file --taxon_output_file phycorder_taxon_stats-$COUNTER-.csv --align_output_file phycorder_total_align_stats-$COUNTER-.csv

    echo $COUNTER > $TEMPFILE
  done

  unlink $TEMPFILE

  # produce data outputs for phycorder alignments

  TEMPFILE=/tmp/$$.tmp
  echo 1 > $TEMPFILE
  for seq_file in $(ls *.aln); do
    COUNTER=$[$(cat $TEMPFILE) + 1]

    $phycorder_path/phy_stats.py --align_file $seq_file --taxon_output_file gon_phy_taxon_stats-$COUNTER-.csv --align_output_file gon_phy_total_align_stats-$COUNTER-.csv

    echo $COUNTER > $TEMPFILE
  done

  unlink $TEMPFILE

  cd $maind

  # move all the valuable files into a folder that shows which run its a part of

  mv ./combined_outputs ./combined_output-$REP_COUNTER
  mv ./phycorder_runs_out ./phycorder_runs_out-$REP_COUNTER
  mv ./phycorder_results ./phycorder_results-$REP_COUNTER
  mv ./gon_phy_runs_dir ./gon_phy_runs_dir-$REP_COUNTER
  mv ./gon_phy_results ./gon_phy_results-$REP_COUNTER

  echo $REP_COUNTER > $REP_TEMPFILE

done

unlink $REP_TEMPFILE

printf "Stepwise phylogenetic analyses finished"
