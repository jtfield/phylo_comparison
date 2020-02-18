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

#TODO: return name of reference for the snippy output

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
		printf "\n$j\n$i\n"
		msa_file=$(basename $j)
		$rapup_path/ref_producer.py -s --align_file $j --out_file $outdir/$gon_phy_basecall/sep_loci/individual_tax/single_tax_$msa_file_$i --ref_select $i
	done
done

