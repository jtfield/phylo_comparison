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
$rapup_path/gon_phyling.sh -d $gon_phy_pwd -1 $r1_tail -2 $r2_tail -o LOCUS -l $loci_positions

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
$rapup_path/multi_map.sh -a $outdir/$update_dir/alignment_ref.fas -d $gon_phy_pwd -1 $r1_tail -2 $r2_tail -o $outdir/rapup_run -g LOCUS -f rapup_loci_positions.csv

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

## SNIPPY FILES BLAST
for j in $(ls $outdir/$loci_blast_indexes/*.fasta | sed -e 's/.fasta//g');
do
	file_name=$(basename $j)
	#printf "\n$file_name\n"
	for i in $(ls $outdir/$snippy_basecall/sep_loci/individual_tax);
	do
		cluster=${i:19:8}
		cluster_grep=$( echo $i | grep -Eo 'cluster[0-9]+')
		#printf "\n$cluster_grep\n"
		#echo ${i:28:14}
		seq=${i:28:14}
		if [ "$seq" == "$file_name" ];
		then
			printf "\nFOUND A MATCH BETWEEN A SNIPPY SEQUENCE AND BLAST INDEX\n"
			blastn -db $j.fasta -query $outdir/$snippy_basecall/sep_loci/individual_tax/$i -out $outdir/$snippy_basecall/blast_results/blast_output_$cluster_grep-$file_name-$seq.out -outfmt 5
			#printf "\n$j\n"
			#printf "\n$i\n"
		fi
	done
done

## RAPUP FILES BLAST
for j in $(ls $outdir/$loci_blast_indexes/*.fasta | sed -e 's/.fasta//g');
do
        file_name=$(basename $j)
        #printf "\n$file_name\n"
        for i in $(ls $outdir/$rapup_basecall/sep_loci/individual_tax);
        do
		cluster=${i:19:8}
		cluster_grep=$( echo $i | grep -Eo 'cluster[0-9]+')
                #printf "\n$cluster_grep\n"i
                #echo ${i:27:14}
                seq=${i:27:14}
                if [ "$seq" == "$file_name" ];
                then
                        printf "\nFOUND A MATCH BETWEEN A RAPUP SEQUENCE AND BLAST INDEX\n"
                        blastn -db $j.fasta -query $outdir/$rapup_basecall/sep_loci/individual_tax/$i -out $outdir/$rapup_basecall/blast_results/blast_output_$cluster_grep-$file_name-$seq.out -outfmt 5
                        #printf "\n$j\n"
                        #printf "\n$i\n"
                fi
        done
done

## gon_phy FILES BLAST
for j in $(ls $outdir/$loci_blast_indexes/*.fasta | sed -e 's/.fasta//g');
do
        file_name=$(basename $j)
        #printf "\n$file_name\n"
        for i in $(ls $outdir/$gon_phy_basecall/sep_loci/individual_tax);
        do
		cluster=${i:19:8}
		cluster_grep=$( echo $i | grep -Eo 'cluster[0-9]+')
                #printf "\n$cluster_grep\n"
                #echo ${i:29:14}
                seq=${i:29:14}
                if [ "$seq" == "$file_name" ];
                then
                        printf "\nFOUND A MATCH BETWEEN A GON_PHY SEQUENCE AND BLAST INDEX\n"
                        blastn -db $j.fasta -query $outdir/$gon_phy_basecall/sep_loci/individual_tax/$i -out $outdir/$gon_phy_basecall/blast_results/blast_output_$cluster_grep-$file_name-$seq.out -outfmt 5
                        #printf "\n$j\n"
                        #printf "\n$i\n"
                fi
        done
done

# 
