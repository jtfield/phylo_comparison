#! /bin/bash 
# pipeline to examine basecall accuracy of RapUp, Snippy and Gon_phyling (de novo assembly) programs

source $1

mkdir -p $outdir

cd $outdir
update_reads=update_reads_dir
start_dir=start_tree_dir
update_dir=update_alignment_dir
loci_blast_indexes=loci_blast_indexes_dir
gon_phy=gon_phyling_dir
blast_output=blast_output_dir
gon_phy_basecall=gon_phy_basecall
rapup_basecall=rapup_basecall
snippy_basecall=snippy_basecall


mkdir $gon_phy

ls "${master_reads}"/*$r1_tail | sort -R > taxa_list.txt

cd $gon_phy

gon_phy_pwd=$(pwd)

for i in $(cat $outdir/taxa_list.txt);
do
	ln -s $i $gon_phy_pwd/
  	ln -s ${i%$r1_tail}$r2_tail $gon_phy_pwd/
done

#CALL GON_PHYLING ON ALL TAXA READS
if [ $intermediate_files == "KEEP" ]; then
	$rapup_path/gon_phyling.sh -d $gon_phy_pwd -1 $r1_tail -2 $r2_tail -o LOCUS -l $loci_positions
elif [ $intermediate_files == "CLEAN" ]; then
	$rapup_path/gon_phyling.sh -d $gon_phy_pwd -1 $r1_tail -2 $r2_tail -o LOCUS -l $loci_positions -i $intermediate_files
fi

sed -i -e '/.ref//g' $gon_phy_pwd/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas
printf "\nGON_PHYLING STAGE COMPLETE\n"
#SET UP REFERENCE BY PULLING A SINGLE SEQUENCE FROM THE GON_PHYLING RUN ALIGNMENT
cd $outdir

mkdir $update_dir

printf "\nBEGINNING REFERENCE SELECTION\n"
$rapup_path/ref_producer.py -r --align_file $gon_phy_pwd/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas --out_file $outdir/$update_dir/alignment_ref.fas

sed -i -e 's/.fasta//g' $outdir/$update_dir/alignment_ref.fas

printf "\nREFERENCE PRODUCTION COMPLETE\n"
ref_name=$(head -1 $outdir/$update_dir/alignment_ref.fas | sed -e 's/>//g' | sed -e 's/.ref//g')

printf "\nref = $ref_name\n"

rm "$gon_phy_pwd/$ref_name"*

ls $gon_phy_pwd

mv $gon_phy_pwd/trimmed_reads $outdir/trimmed_reads

#BEGIN ROUNDS OF RAPUP AND SNIPPY
if [ $intermediate_files == "KEEP" ]; then
	$rapup_path/multi_map.sh -a $outdir/$update_dir/alignment_ref.fas -d $gon_phy_pwd -1 $r1_tail -2 $r2_tail -o $outdir/rapup_run -g LOCUS -f rapup_loci_positions.csv
elif [ $intermediate_files == "CLEAN" ]; then
	$rapup_path/multi_map.sh -a $outdir/$update_dir/alignment_ref.fas -d $gon_phy_pwd -1 $r1_tail -2 $r2_tail -o $outdir/rapup_run -g LOCUS -f rapup_loci_positions.csv -i $intermediate_files
fi

touch $update_dir/snippy_run.tab

#SET UP FOR SNIPPY MULTI-TAXON RUN
for i in $(ls $gon_phy_pwd/*$r1_tail);
do
	isolate_name=$(basename $i $r1_tail)
	printf "$isolate_name"
	printf "\n$i\n"
	printf "\n${i%$r1_tail}$r2_tail\n"
	echo "$isolate_name	$i	${i%$r1_tail}$r2_tail" >> $update_dir/snippy_run.tab

done

$snippy_path/snippy-multi $outdir/$update_dir/snippy_run.tab --ref $outdir/$update_dir/alignment_ref.fas --cpus 4 > $outdir/runme.sh

chmod +x $outdir/runme.sh

$outdir/runme.sh

# WARNING WARNING WARNING: THIS RAXML COMMAND IS WONKY AS FUCK
#raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $outdir/core.full.aln -p 12345 -n consensusFULL

#This one works for some ungodly reason. wtf and why cant i call -T for threads?
raxmlHPC-PTHREADS -s $outdir/core.full.aln -m GTRGAMMA -p 12345 -n snippy_tree


#trying to clean up snippy...
if [ $intermediate_files == "CLEAN" ]; then
        printf "\nCleaning up intermediate output files.\n"
        for i in $(cat $outdir/taxa_list.txt); do
		base_file=$(basename $i $r1_tail)
		printf "\n$base_file\n"
		rm -r $base_file
                #cd $i
                #for j in $(ls -1); do
                #        rm ./$j
                #done    
                #cd ..
                #rmdir $i
        done
elif [ $intermediate_files == "KEEP" ]; then
        printf "\nKeeping intermediate output files\n"
fi



#BEGIN SEPARATION OF LOCI IF NECESSARY AND ALIGNMENT OF CONSTRUCTED LOCI TO THE ORIGINAL GENOMES THAT PRODUCED THE READS

# CONSTRUCT BLAST INDEXES
mkdir $outdir/$loci_blast_indexes

for i in $(ls $reference_dir);
do
	printf "\n$reference_dir/$i\n"
	ln -s $reference_dir/$i $outdir/$loci_blast_indexes/$i
	makeblastdb -in $outdir/$loci_blast_indexes/$i -dbtype nucl -parse_seqids

done

mkdir $outdir/$rapup_basecall
mkdir $outdir/$snippy_basecall
mkdir $outdir/$gon_phy_basecall
mkdir $outdir/$rapup_basecall/blast_results
mkdir $outdir/$snippy_basecall/blast_results
mkdir $outdir/$gon_phy_basecall/blast_results

#SPLIT EACH LOCUS INTO A SEPARATE ALIGNMENT FILE

#TODO: thuroughly check replacement of name of reference for the snippy output

#REPLACE THE term REFERENCE with the actual reference name in the snippy output alignment
sed -i -e 's/Reference/'$ref_name'/g' $outdir/core.full.aln

$rapup_path/locus_position_identifier.py --out_file_dir $outdir/$rapup_basecall/sep_loci --position_csv_file $loci_positions --concatenated_fasta $outdir/rapup_run/combine_and_infer/extended.aln

$rapup_path/locus_position_identifier.py --out_file_dir $outdir/$snippy_basecall/sep_loci --position_csv_file $loci_positions --concatenated_fasta $outdir/core.full.aln

$rapup_path/locus_position_identifier.py --out_file_dir $outdir/$gon_phy_basecall/sep_loci --position_csv_file $loci_positions --concatenated_fasta $outdir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas

mkdir $outdir/$rapup_basecall/sep_loci/individual_tax
mkdir $outdir/$snippy_basecall/sep_loci/individual_tax
mkdir $outdir/$gon_phy_basecall/sep_loci/individual_tax

for j in $(ls $outdir/$gon_phy_basecall/sep_loci/*.fasta);
do
	for i in $(grep ">" $j | sed -e 's/>//g');
	do
		#printf "\n$j\n$i\n"
		msa_file=$(basename $j | sed -e 's/.fasta//g')
		#printf "\n$msa_file\n"
		$rapup_path/ref_producer.py -s --align_file $j --out_file $outdir/$gon_phy_basecall/sep_loci/individual_tax/single_tax_gon_phy-$msa_file-$i --ref_select $i
	done
done

for j in $(ls $outdir/$rapup_basecall/sep_loci/*.fasta);
do
        for i in $(grep ">" $j | sed -e 's/>//g');
        do
                #printf "\n$j\n$i\n"
                msa_file=$(basename $j | sed -e 's/.fasta//g')
		#printf "\n$msa_file\n"
                $rapup_path/ref_producer.py -s --align_file $j --out_file $outdir/$rapup_basecall/sep_loci/individual_tax/single_tax_rapup-$msa_file-$i --ref_select $i
        done
done

#Create single taxa locus files for snippy output
for j in $(ls $outdir/$snippy_basecall/sep_loci/*.fasta);
do
        for i in $(grep ">" $j | sed -e 's/>//g');
        do
                #printf "\n$j\n$i\n"
                msa_file=$(basename $j | sed -e 's/.fasta//g')
		#printf "\n$msa_file\n"
                $rapup_path/ref_producer.py -s --align_file $j --out_file $outdir/$snippy_basecall/sep_loci/individual_tax/single_tax_snippy-$msa_file-$i --ref_select $i
        done
done

printf "\n BEGINNING BLAST PROCESS\n"


#FIND THE CORRECT SEQUENCES TO BLAST AGAINST THE TRUE GENOME SEQUENCES
#TODO: the semi-hardcoding of file name based on structure isnt good coding and needs to be changed

touch $outdir/$snippy_basecall/snippy_alignment_runs.txt

cd $outdir/$snippy_basecall/sep_loci/individual_tax

## SNIPPY FILES BLAST
for j in $(ls $outdir/$loci_blast_indexes/*.fasta | sed -e 's/.fasta//g');
do
	file_name=$(basename $j)
	printf "\n$file_name\n"
	for i in $(ls $outdir/$snippy_basecall/sep_loci/individual_tax);
	do
		#cluster=${i:19:8}
		cluster_grep=$( echo $i | grep -Eo 'cluster[0-9]+')
		#printf "\n$cluster_grep\n"
		#echo ${i:28:14}
		seq_grep=$( echo $i | grep -o "$file_name$")
		#printf "\n $seq_grep \n"
		#printf '\nFILE TO BE MATCHED == '$i'\n'
		#seq=${i:28:14}
		#if [ "$seq" == "$file_name" ];
		if [ ! -z "$seq_grep" ];
		then
			printf "\nFOUND A MATCH BETWEEN A SNIPPY SEQUENCE AND BLAST INDEX\n"
			
			$program_path/seq_producer.py --align $outdir/$snippy_basecall/sep_loci/individual_tax/$i --output_align_stub $cluster_grep-$seq_grep

			cat $cluster_grep-$seq_grep-complement.fasta <(echo) $j.fasta > combine-$cluster_grep-$seq_grep-complement.fasta
			cat $cluster_grep-$seq_grep-reverse_complement.fasta <(echo) $j.fasta > combine-$cluster_grep-$seq_grep-reverse_complement.fasta
			cat $cluster_grep-$seq_grep-reverse.fasta <(echo) $j.fasta > combine-$cluster_grep-$seq_grep-reverse.fasta
			cat $i <(echo) $j.fasta > combine-$cluster_grep-$seq_grep-original.fasta
			
			mafft combine-$cluster_grep-$seq_grep-complement.fasta > aligned_combine-$cluster_grep-$seq_grep-complement.fasta
			mafft combine-$cluster_grep-$seq_grep-reverse_complement.fasta > aligned_combine-$cluster_grep-$seq_grep-reverse_complement.fasta
			mafft combine-$cluster_grep-$seq_grep-reverse.fasta > aligned_combine-$cluster_grep-$seq_grep-reverse.fasta
			mafft combine-$cluster_grep-$seq_grep-original.fasta > aligned_combine-$cluster_grep-$seq_grep-original.fasta


			#Perform basecall comparison on all method sequences with the sequence that produced the reads
			$program_path/align_compare.py --align_1 aligned_combine-$cluster_grep-$seq_grep-reverse_complement.fasta --align_2 aligned_combine-$cluster_grep-$seq_grep-complement.fasta --align_3 aligned_combine-$cluster_grep-$seq_grep-reverse.fasta --align_4 aligned_combine-$cluster_grep-$seq_grep-original.fasta --output_stub basecall_results-$cluster_grep-$seq_grep-.txt

			#echo "$program_path/alignment_checker.py --align_1 $j.fasta --align_2 $outdir/$snippy_basecall/sep_loci/individual_tax/$i --output_align $outdir/$snippy_basecall/blast_results/realignment_$cluster_grep-$file_name-$seq_grep.fasta --output_miscalls $outdir/$snippy_basecall/blast_results/miscalls_output_$cluster_grep-$file_name-$seq_grep.out" >> $outdir/$snippy_basecall/snippy_alignment_runs.txt
			
			#echo "$program_path/alignment_checker.py+--align_1+$j.fasta+--align_2+$outdir/$snippy_basecall/sep_loci/individual_tax/$i+--output_align+$outdir/$snippy_basecall/blast_results/realignment_$cluster_grep-$file_name-$seq_grep.fasta+--output_miscalls+$outdir/$snippy_basecall/blast_results/miscalls_output_$cluster_grep-$file_name-$seq_grep.out" >> $outdir/$snippy_basecall/snippy_alignment_runs.txt


			#blastn -db $j.fasta -query $outdir/$snippy_basecall/sep_loci/individual_tax/$i -out $outdir/$snippy_basecall/blast_results/blast_output_$cluster_grep-$file_name-$seq_grep.out -max_hsps 1 -outfmt 5
			#printf "\n$j\n"
			#printf "\n$i\n"
		fi
	done
done

cd $outdir

#split up alignment runs for parallel processing
#cd $outdir/$snippy_basecall/

#cat $outdir/$snippy_basecall/snippy_alignment_runs.txt | split -d -l $align_threads

#for j in $(ls x*); do
#        for i in $(cat $j); do
#        	
#		echo "$i"
#                echo "split"
#                split=$(echo $i | sed -e 's/+/ /g')
#                echo "$split"
#                echo "waffle"
#		split_array=($split)
#		$program_path/parallel_alignment_checker.py --align_1 ${split_array[2]} --align_2 ${split_array[4]} --output_align ${split_array[6]} --output_miscalls ${split_array[8]}
#		echo "${split_array[0]}"
#                echo "waffle"
#
#
#	done
#        wait
#done
#
#cd -



## RAPUP FILES BLAST

touch $outdir/$rapup_basecall/rapup_alignment_runs.txt

cd $outdir/$rapup_basecall/sep_loci/individual_tax

for j in $(ls $outdir/$loci_blast_indexes/*.fasta | sed -e 's/.fasta//g');
do
        file_name=$(basename $j)
        printf "\n$file_name\n"
        for i in $(ls $outdir/$rapup_basecall/sep_loci/individual_tax);
        do
		#cluster=${i:19:8}
		cluster_grep=$( echo $i | grep -Eo 'cluster[0-9]+')
		seq_grep=$( echo $i | grep -o "$file_name$")
                #printf "\n$cluster_grep\n"i
                #echo ${i:27:14}
                #seq=${i:27:14}
                #if [ "$seq" == "$file_name" ];
		if [ ! -z "$seq_grep" ];
                then
                        printf "\nFOUND A MATCH BETWEEN A RAPUP SEQUENCE AND BLAST INDEX\n"
                        $program_path/seq_producer.py --align $outdir/$rapup_basecall/sep_loci/individual_tax/$i --output_align_stub $cluster_grep-$seq_grep
			
			#echo "$program_path/alignment_checker.py --align_1 $j.fasta --align_2 $outdir/$rapup_basecall/sep_loci/individual_tax/$i --output_align $outdir/$rapup_basecall/blast_results/realignment_$cluster_grep-$file_name-$seq_grep.fasta --output_miscalls $outdir/$rapup_basecall/blast_results/miscalls_output_$cluster_grep-$file_name-$seq_grep.out" >> $outdir/$rapup_basecall/rapup_alignment_runs.txt
			
			#blastn -db $j.fasta -query $outdir/$rapup_basecall/sep_loci/individual_tax/$i -out $outdir/$rapup_basecall/blast_results/blast_output_$cluster_grep-$file_name-$seq_grep.out -max_hsps 1 -outfmt 5
                        #printf "\n$j\n"
                        #printf "\n$i\n"
                fi
        done
done

cd $outdir

#split up alignment runs for parallel processing
#cd $outdir/$rapup_basecall/

#cat $outdir/$rapup_basecall/rapup_alignment_runs.txt | split -d -l $align_threads

#for j in $(ls x*); do
#        #for i in $(cat $j  | ls -1); do
#        while read i; do
#               echo "${i[0]}"
#        done < $j
#        wait
#done

#cd -



## gon_phy FILES BLAST

cd $outdir/$gon_phy_basecall/sep_loci/individual_tax

touch $outdir/$gon_phy_basecall/gon_phy_alignment_runs.txt

for j in $(ls $outdir/$loci_blast_indexes/*.fasta | sed -e 's/.fasta//g');
do
        file_name=$(basename $j)
        printf "\n$file_name\n"
        for i in $(ls $outdir/$gon_phy_basecall/sep_loci/individual_tax);
        do
		#cluster=${i:19:8}
		cluster_grep=$( echo $i | grep -Eo 'cluster[0-9]+')
                seq_grep=$( echo $i | grep -o "$file_name$")
		#printf "\n$cluster_grep\n"
                #echo ${i:29:14}
                #seq=${i:29:14}
                #if [ "$seq" == "$file_name" ];
		if [ ! -z "$seq_grep" ];
                then
                        printf "\nFOUND A MATCH BETWEEN A GON_PHY SEQUENCE AND BLAST INDEX\n"
                        $program_path/seq_producer.py --align $outdir/$gon_phy_basecall/sep_loci/individual_tax/$i --output_align_stub $cluster_grep-$seq_grep
			
			#echo "$program_path/alignment_checker.py --align_1 $j.fasta --align_2 $outdir/$gon_phy_basecall/sep_loci/individual_tax/$i --output_align $outdir/$gon_phy_basecall/blast_results/realignment_$cluster_grep-$file_name-$seq_grep.fasta --output_miscalls $outdir/$gon_phy_basecall/blast_results/miscalls_output_$cluster_grep-$file_name-$seq_grep.out" >> $outdir/$gon_phy_basecall/gon_phy_alignment_runs.txt
			
			#blastn -db $j.fasta -query $outdir/$gon_phy_basecall/sep_loci/individual_tax/$i -out $outdir/$gon_phy_basecall/blast_results/blast_output_$cluster_grep-$file_name-$seq_grep.out -max_hsps 1 -outfmt 5
                        #printf "\n$j\n"
                        #printf "\n$i\n"
                fi
        done
done

cd $outdir

#split up alignment runs for parallel processing
#cd $outdir/$gon_phy_basecall/

#cat $outdir/$gon_phy_basecall/gon_phy_alignment_runs.txt | split -d -l $align_threads

#for j in $(ls x*); do
#        #for i in $(cat $j  | ls -1); do
#	while read i; do
#               echo "${i[0]}"
#        done < $j
#        wait
#done

#cd -
# COPY TRUE TREE TO OUTPUT FOLDER
cp $true_tree $outdir/true_tree.tre

printf "\n$program_path\n"
