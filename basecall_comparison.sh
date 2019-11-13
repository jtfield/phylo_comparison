#!/bin/bash
# Phycorder and Gon_phyling basecall accuracy comparison
PHY_COMPARE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


source $1

mkdir -p $outdir

cd $outdir

# split read names for the phycorder starting tree into a file
ls "${master_reads}"/*$r1_tail | head -2 > start_tree_files.txt

# split reads for the rest of the updating using phycorder into a seperate file
ls "${master_reads}"/*$r1_tail | tail -n +3 > update_read_list.txt

start_dir=start_tree_dir
update_dir=update_alignment_dir
loci_blast_indexes=loci_blast_indexes_dir
truncated_tree=phycorder_low_loci_dir
gon_phy=gon_phyling_dir
gon_phy_tree=gon_phyling_tree_dir
blast_output=blast_output_dir
gon_phy_basecall=gon_phy_basecall
phycord_basecall=phycord_basecall


mkdir $start_dir
mkdir $truncated_tree
mkdir $update_dir
mkdir $loci_blast_indexes
mkdir $gon_phy
mkdir $gon_phy_tree
mkdir $blast_output
mkdir $gon_phy_basecall
mkdir $phycord_basecall


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

for i in $(cat $outdir"/"start_tree_files.txt); do
  ln -s "${read_dir}"/$i $gon_phy
  ln -s "${read_dir}"/${i%$r1_tail}$r2_tail $gon_phy
done

# link files for updating the starting alignment to a seperate directory
for i in $(cat $outdir"/"update_read_list.txt); do
  ln -s "${read_dir}"/$i $gon_phy
  ln -s "${read_dir}"/${i%$r1_tail}$r2_tail $gon_phy
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

printf "USING THIS NUMBER OF LOCI == $num_loci"

$PHY_COMPARE/select_ten.py --msa_folder $outdir/$start_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/locus_msa_files --position_csv_file $loci_positions --out_file $outdir/$loci_blast_indexes/loci_for_use.txt --num_loci $num_loci

# move a set number (max of 10) of loci to a folder where blast indexes will be constructed
# SEPERATE TOP TAXON SEQUNCE AS REPRESENTITIVE TO BE USED FOR BLAST INDEX
for i in $(cat $outdir"/"$loci_blast_indexes"/"loci_for_use.txt); do

	cp "$outdir/$start_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/locus_msa_files/$i" $outdir/$loci_blast_indexes/

	#cp "$outdir/$start_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/locus_msa_files/$i" $outdir/$truncated_tree/

	$phycorder_path/ref_producer.py -s --align_file $outdir/$loci_blast_indexes/$i --out_file $outdir/$loci_blast_indexes/$i-single.fasta

	makeblastdb -in $outdir/$loci_blast_indexes/$i-single.fasta -dbtype nucl -parse_seqids

done

# COMBINE SELECTED LOCI INTO SINGLE FASTA FILE FOR TRUNCATED PHYCORDER RUN
#$phycorder_path/locus_combiner.py --msa_folder $outdir/$truncated_tree/ --out_file $outdir/$truncated_tree/phycord_base.fasta --position_dict_file $outdir/$truncated_tree/pos_file.txt --suffix .fasta

# RUN GON_PHYLING ON READS FOR ALL TAXA IN LOCUS MODE
if [ $gon_phy_alignments == "OFF" ]; then

# create a config file on the fly with current settings
  cat <<EOF > gon_phy_assembly.cfg
  ## Config file for gon_phyling pipeline
  # change the values of the variables to control the pipeline
  # the variables are structured like this:
  # variable_name="value_of_the_variable"
  # value of the variable can be a path to a file or directory of files
  # Path to the reference genome required for Parsnp
  ref_genome="$ref_genome"
  # Path to the directory of reads for assembly
  read_dir="$outdir/$gon_phy"
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
  #loci_positions="$loci_positions"
  loci_positions="$outdir/gon_phy_run_positions.txt"
EOF

fi

time $phycorder_path/gon_phyling.sh $outdir/gon_phy_assembly.cfg

for j in $(ls $outdir/$loci_blast_indexes/*single.fasta); do
        for i in $(ls $outdir/$gon_phy/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/locus_msa_files/*.fasta); do

                database=$(basename $j)
                query=$(basename $i)
                blastn -db $j -query $i -out $blast_output/blast_output-$database-$query.txt -outfmt 5
                #printf "$j \n"
                #printf "$i \n"
                #printf "\n"



        done
        wait
done


$PHY_COMPARE/basecall_loci_matcher.py --blast_output_folder $outdir/$blast_output --gon_out_file $outdir/gon_phyling_files.txt --phy_out_file $outdir/phycorder_files.txt


# MOVE FILES INTO APPROPRIATE FOLDERS

for i in $(cat $outdir"/"phycorder_files.txt); do
	cp "$outdir/$start_dir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/locus_msa_files/$i" $outdir/$truncated_tree/
	printf "CLUSTER FOR PHYCORDER $i"
done

#$phycorder_path/locus_combiner.py --msa_folder $outdir/$truncated_tree/ --out_file $outdir/$truncated_tree/phycord_base.fasta --position_dict_file $outdir/$truncated_tree/phycord_pos_file.txt --suffix .fasta

echo "phycorder_path/new_locus_combiner.py --msa_folder $outdir/$truncated_tree/ --suffix .fasta --out_file $outdir/$truncated_tree/phycord_base.fasta --len_filter 1000 --position_csv_file $outdir/$truncated_tree/phycord_pos_file.csv"


$phycorder_path/new_locus_combiner.py --msa_folder $outdir/$truncated_tree/ --suffix .fasta --out_file $outdir/$truncated_tree/phycord_base.fasta --len_filter 1000 --position_csv_file $outdir/$truncated_tree/phycord_pos_file.csv


#for i in $(ls $outdir/$gon_phy/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/locus_msa_files/*.fasta); do
for i in $(cat $outdir"/"gon_phyling_files.txt); do
	 cp $outdir/$gon_phy/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/locus_msa_files/$i /$outdir/$gon_phy_tree/

done

#$phycorder_path/locus_combiner.py --msa_folder $outdir/$gon_phy_tree/ --out_file $outdir/$gon_phy_tree/gon_phy_base.fasta --position_dict_file $outdir/$gon_phy_tree/gon_phypos_file.txt --suffix .fasta

echo "$phycorder_path/new_locus_combiner.py --msa_folder $outdir/$gon_phy_tree/ --suffix .fasta --out_file $outdir/$gon_phy_tree/gon_phy_base.fasta --len_filter 1000 --position_csv_file $outdir/$gon_phy_tree/gon_phypos_file.csv"

$phycorder_path/new_locus_combiner.py --msa_folder $outdir/$gon_phy_tree/ --suffix .fasta --out_file $outdir/$gon_phy_tree/gon_phy_base.fasta --len_filter 1000 --position_csv_file $outdir/$gon_phy_tree/gon_phypos_file.csv



# RUN RAxML ON GON_PHYLING LOCI
cd $outdir/$gon_phy_tree

time raxmlHPC-PTHREADS -m GTRGAMMA -T $THREADS -s $outdir/$gon_phy_tree/gon_phy_base.fasta -p 12345 -n core_genome_run.out

cd $outdir
# RUN PHYCORDER ON STARTING TREE LOCI

# MAKE CONFIG FILE FOR PHYCORDER
#cat <<phy_loop > basic.cfg

        ## Welcome to the Phycorder config file
        # Change the variable values to match the files and numbers you wish to use
        # this is the path to the phycorder directory. Dont move it.

        # PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

        # the path to the previously generated alignment file

#        align="$outdir/$truncated_tree/phycord_base.fasta"

        # the path to the previously generated tree file made from the alignment file

#        tree="NONE"

        # the directory of paired end read pairs you belonging to taxa you wish to add to your phylogeny

#        read_dir="$outdir/$update_dir"

        # the number of taxa that can be added to the phylogeny at a single time
        # should be less than the number of cores available
        # should be balanced with the number of threads you will assign to the programs within phycorder

#        phycorder_runs="$phycorder_runs"

        # the number of threads you wish to make available to each phycorder run
        # for mapping with bowtie and inference with RAxML

#        threads="$phycorder_threads"

        # the tail identifiers of the read pairs
        # if the full read name is "Clade_1_01_R1_001.fastq" and Clade_1_01_R2_001.fastq"
        # then only put the portion of the file names that change to correspond to the read pairs
        # in this example, Clade_1_01_ identify the taxons and so must not be included

#        r1_tail="$r1_tail"
#        r2_tail="$r2_tail"

        # the output directory for your final information

#        outdir="$outdir/phycorder-out"

        #bootstrapping
#        bootstrapping=$bootstrapping

        # In addition, you'll currently need to run Phycorder with align_type="LOCUS" to run this mode
        # you'll need to pass in the location of the output positional dictionary file describing the length
        # of the loci in the concatenated alignment. This is important because it also contains the relative positon
        # in the concatenated alignment
        # give this variable the absolute path to the dict file
#        loci_positions="$outdir/$truncated_tree/phycord_pos_file.txt"

        # depending on how you want to use phycorder, you'll change this variable value.
        # if you have multiple sequence alignments for multiple single locus,
        # if you have a fasta file with concatenated MSA and a seperate file containing loci start and stop positions
        # if you have a nexus file with concatenated MSA and loci start and stop positions
        # EXAMPLE: align_type="NEXUS" NOT SUPPORTED YET<<<<<<<<<<<<<<<<<<
        # EXAMPLE: align_type="SINGLE_LOCUS_FILES"
        # EXAMPLE: align_type="CONCAT_MSA"
        # EXAMPLE: align_type="PARSNP_XMFA"
#        align_type="CONCAT_MSA"
        #align_type="PARSNP_XMFA"
        #align_type="SINGLE_LOCUS_FILES"

	# depending on how you want to use phycorder, you'll change this variable value.
        # if you want multiple sequence alignments for multiple single locus,
        # if you want  a fasta file with concatenated MSA and a seperate file containing loci start and stop positions
        # if you want a nexus file with concatenated MSA and loci start and stop positions
        # EXAMPLE: output_type="NEXUS" NOT SUPPORTED YET<<<<<<<<<<<<<<<<<<<<
        # EXAMPLE: output_type="SINGLE_LOCUS_FILES"
        # EXAMPLE: output_type="CONCAT_MSA"
#        output_type="SINGLE_LOCUS_FILES"
        #output_type="CONCAT_MSA"

#phy_loop

#$phycorder_path/multi_map.sh $outdir/basic.cfg


# commandline phycorder

printf $phycorder_path/multi_map.sh -a "$outdir/$truncated_tree/phycord_base.fasta" -d "$outdir/$update_dir" -m CONCAT_MSA -g SINGLE_LOCUS_FILES -1 "$r1_tail" -2 "$r2_tail" -p $phycorder_runs -c $phycorder_threads -o "$outdir/phycorder-out" -f "$outdir/$truncated_tree/phycord_pos_file.csv"



$phycorder_path/multi_map.sh -a "$outdir/$truncated_tree/phycord_base.fasta" -d "$outdir/$update_dir" -m CONCAT_MSA -g SINGLE_LOCUS_FILES -1 "$r1_tail" -2 "$r2_tail" -p $phycorder_runs -c $phycorder_threads -o "$outdir/phycorder-out" -f "$outdir/$truncated_tree/phycord_pos_file.csv"


# MOVE ALL OUTPUT LOCUS FILES TO THE APPROPRIATE METHOD BASECALL FOLDERS
ln -s $outdir/phycorder-out/combine_and_infer/updated_single_loci/*.fasta $outdir/phycord_basecall

ln -s $outdir/gon_phyling_tree_dir/cluster*.fasta $outdir/gon_phy_basecall

# SEPARATE OUT EACH TAXON SEQUENCE FOR EACH LOCUS TO ITS OWN FILE

for i in $(ls -l $outdir/phycord_basecall); do

	$PHY_COMPARE/sequence_separate.py --concatenated_fasta $i --out_file_dir $outdir/phycord_basecall

done

for i in $(ls -l $outdir/gon_phy_basecall); do

        $PHY_COMPARE/sequence_separate.py --concatenated_fasta $i --out_file_dir $outdir/gon_phy_basecall

done

# BLAST ALL SEPARATED TAXON LOCI AGAINST THE GENOMES THAT PRODUCED THEM

printf "FILES READY FOR BLAST INDEX CREATION"

mkdir final_blast_idxs

for i in $(ls $reference_dir); do

	ln -s $reference_dir/$i $outdir/final_blast_idxs/
	makeblastdb -in $outdir/final_blast_idxs/$i -dbtype nucl -parse_seqids

done

# CALL PROGRAM TO MATCH SEQUENCE FILE SETS FROM BOTH METHODS TO APPROPRIATE REFERENCE FILES
# OUTPUT FILES WITH THE SEQUENCE QUERY FILE NAME AND THE REFERENCE FILE NAME TO RESPECTIVE FOLDERS

$PHY_COMPARE/basecall_match_blast.py --ref_seq_dir $outdir/final_blast_idxs --seq_set_1_dir $outdir/gon_phy_basecall --seq_set_2_dir $outdir/phycord_basecall --ref_suffix .fasta

# MAKE FOLDERS FOR BLAST RESULTS OF BOTH METHODS

mkdir $outdir/final_phycord_blast_output
mkdir $outdir/final_gon_phy_blast_output

# RUN BLAST AND OUTPUT TO RESPECTIVE FOLDERS
for i in $(ls $outdir/gon_phy_basecall/*.txt); do
	base=$(basename $i)
	query=$(head -1 $i)
	ref=$(tail -1 $i)
	#blastn -db $ref -query $query -out $outdir/final_gon_phy_blast_output/$base-blast-out.txt -outfmt 5
	blastn -db $ref -query $query -out $outdir/final_gon_phy_blast_output/$base-blast-out.txt -outfmt "10 qstart qend sstart send length pident nident mismatch gapopen gaps"  

done

# RUN BLAST AND OUTPUT TO RESPECTIVE FOLDERS
for i in $(ls $outdir/phycord_basecall/*.txt); do
        base=$(basename $i)
        query=$(head -1 $i)
        ref=$(tail -1 $i)
        #blastn -db $ref -query $query -out $outdir/final_phycord_blast_output/$base-blast-out.txt -outfmt 5
	blastn -db $ref -query $query -out $outdir/final_phycord_blast_output/$base-blast-out.txt -outfmt "10 qstart qend sstart send length pident nident mismatch gapopen gaps"


done

printf "BLASTING OF METHOD CONSTRUCTED SEQUENCES AGAINST SEQUENCES THEY WERE SIMULATED FROM IS COMPLETE"
