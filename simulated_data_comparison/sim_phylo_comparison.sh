#! /bin/bash 
# pipeline to examine basecall accuracy of RapUp, Snippy and Gon_phyling (de novo assembly) programs

source $1
if [[ $assembly == "ON" ]];
then 

	update_reads=update_reads_dir
	start_dir=start_tree_dir
	update_dir=update_alignment_dir
	loci_blast_indexes=loci_blast_indexes_dir
	gon_phy=gon_phyling_dir
	blast_output=blast_output_dir
	gon_phy_basecall=gon_phy_basecall
	rapup_basecall=rapup_basecall
	snippy_basecall=snippy_basecall
	snippy_rapup_basecall=snip_rap_basecall
	short_ref=five_align
	snippy_reads=snippy_reads

	mkdir -p $outdir

	printf "\n$outdir\n"

	cd $outdir

	mkdir $snippy_reads
	mkdir $gon_phy

	cd $gon_phy

	gon_phy_pwd=$(pwd)

	cd $outdir

	cd $snippy_reads
	snippy_reads_pwd=$(pwd)

	cd ${gon_phy_pwd}


	# Make folders holding reads for gon_phyling/rapup and Snippy
	# This folder is for gon_phyling and rapup
	ls ${master_reads}/*$r1_tail | sort -R > ${outdir}/taxa_list.txt

	for i in $(cat ${outdir}/taxa_list.txt);
	do
		ln -s $i $gon_phy_pwd/
	  	ln -s ${i%$r1_tail}$r2_tail $gon_phy_pwd/
	done

	# this folder is for Snippy
	for i in $(cat ${outdir}/taxa_list.txt);
	do
		ln -s $i ${snippy_reads_pwd}/
	  	ln -s ${i%$r1_tail}$r2_tail ${snippy_reads_pwd}/
	done

	#CALL GON_PHYLING ON ALL TAXA READS
	if [ $intermediate_files == "KEEP" ]; then
		${rapup_path}/gon_phyling.sh -d $gon_phy_pwd -1 $r1_tail -2 $r2_tail -o LOCUS -l $loci_positions -b $bootstrapping
	elif [ $intermediate_files == "CLEAN" ]; then
		${rapup_path}/gon_phyling.sh -d $gon_phy_pwd -1 $r1_tail -2 $r2_tail -o LOCUS -l $loci_positions -i $intermediate_files -b $bootstrapping
	fi

	sed -i "s/.ref//g" ${gon_phy_pwd}/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas

	printf "\nGON_PHYLING STAGE COMPLETE\n"

	# Run parsnp on a small subset of sequences to use a sequence as a reference for rapup
	# This is so we don't bias our rapup run with baked-in knowledge of indels and missing data
	# from a complete gon_phyling run

	printf "\nBeginning small parsnp run for rapup reference.\n"
	mkdir ${outdir}/${short_ref}

	snip_ref_file=$(echo "${snip_ref}.fasta")
	cp ${gon_phy_pwd}/trimmed_reads/spades_output/genomes_for_parsnp/${snip_ref_file} ${outdir}/${short_ref}/

	rm ${gon_phy_pwd}/${snip_ref}${r1_tail}
	rm ${gon_phy_pwd}/${snip_ref}${r2_tail}

	for i in $(ls -1 ${gon_phy_pwd}/trimmed_reads/spades_output/genomes_for_parsnp/*.fasta | sort -R | head -4);
	do
		ref_read_names=$(basename $i .fasta)
		echo "$ref_read_names"
		echo "$i"
		echo "$snip_ref"
		if [[ $ref_read_names != ${snip_ref} ]];
		then

			cp ${i} ${outdir}/${short_ref}/
			rm ${gon_phy_pwd}/${ref_read_names}${r1_tail}
			rm ${gon_phy_pwd}/${ref_read_names}${r2_tail}
			printf "\n${i}\n"
		elif [[ $ref_read_names == ${snip_ref} ]];
		then
			echo "Skipping reference assembly as this should already be in the five_align folder"
			# replace_ref=$(ls -1 ${gon_phy_pwd}/trimmed_reads/spades_output/genomes_for_parsnp/*.fasta | tail -2 | head -1)
			# ref_read_name=$(basename $replace_ref .fasta)
			# cp ${replace_ref} ${outdir}/${short_ref}/
			# rm ${gon_phy_pwd}/${ref_read_names}${r1_tail}
			# rm ${gon_phy_pwd}/${ref_read_names}${r2_tail}
			# printf "\n${replace_ref}\n"

		fi
	done

	cd ${outdir}/${short_ref}

	parsnp -c -d ${outdir}/${short_ref} -r !

	sed -i 's/.fasta//g' ${outdir}/${short_ref}/P*/parsnp.xmfa
	sed -i 's/.ref//g' ${outdir}/${short_ref}/P*/parsnp.xmfa

	# snip_ref=$(ls -1 | head -1 | sed 's/.fasta//g')

	printf "\nref = $snip_ref\n"

	# rm "${gon_phy_pwd}/$snip_ref"*

	${rapup_path}/multi_map.sh \
	-a ${outdir}/${short_ref}/P*/parsnp.xmfa \
	-r $snip_ref \
	-m PARSNP_XMFA \
	-d $gon_phy_pwd \
	-1 $r1_tail -2 $r2_tail \
	-o ${outdir}/rapup_run \
	-g LOCUS \
	-f rapup_loci_positions.csv \
	-i $intermediate_files \
	-b $bootstrapping

	cd ${outdir}/rapup_run/combine_and_infer

	# sed -i 's/.fasta//g' extended.aln
	# sed -i 's/.ref//g' extended.aln

	#SET UP REFERENCE BY PULLING A SINGLE SEQUENCE FROM THE GON_PHYLING RUN ALIGNMENT
	cd $outdir

	mkdir $update_dir

	####################################################################################
	#SETUP FOR SNIPPY RUN
	printf "\nFINDING REFERENCE FOR SNIPPY AND RAPUP\n"
	selected_ref=$(grep "$snip_ref" ${gon_phy_pwd}/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas | sed -e 's/>//g')
	printf "\n$selected_ref\n"
	#snippy_ref_folder=$(ls $outdir/gon_phyling_dir/trimmed_reads/spades_output/*$snip_ref*)

	#snippy_ref=${outdir}/gon_phyling_dir/trimmed_reads/spades_output/genomes_for_parsnp/$snip_ref.fasta
	snippy_ref=${ref_folder}/${snip_ref}.fasta

	printf "\n$snippy_ref\n"

	cp $snippy_ref $outdir/$update_dir/snippy_ref.fas

	# sed -i -e 's/>.+$/'$snip_ref'/g' $outdir/$update_dir/snippy_ref.fas

	#printf "\nBEGINNING REFERENCE SELECTION\n"
	#$rapup_path/modules/ref_producer.py -s --ref_select $selected_ref --align_file $gon_phy_pwd/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas --out_file $outdir/$update_dir/alignment_ref.fas

	touch ${update_dir}/snippy_run.tab

	#SET UP FOR SNIPPY MULTI-TAXON RUN
	for i in $(ls ${snippy_reads_pwd}/*$r1_tail);
	do
		isolate_name=$(basename $i $r1_tail)
		if [ $isolate_name != $snip_ref ]; then
			printf "$isolate_name"
			printf "\n$i\n"
			printf "\n${i%$r1_tail}$r2_tail\n"
			echo "$isolate_name	$i	${i%$r1_tail}$r2_tail" >> ${update_dir}/snippy_run.tab
		fi

	done

	snippy-multi ${outdir}/${update_dir}/snippy_run.tab --ref ${outdir}/${update_dir}/snippy_ref.fas --cpus 4 > ${outdir}/runme.sh

	chmod +x $outdir/runme.sh

	${outdir}/runme.sh

	#REPLACE THE term REFERENCE with the actual reference name in the snippy output alignment
	printf "\nREPLACE SNIPPY REF\n"
	sed -i -e 's/Reference/'$snip_ref'/g' $outdir/core.full.aln

########################################################################################################################
#SNIPPY PHYLOGENY INFERENCE
	# # WARNING WARNING WARNING: THIS RAXML COMMAND IS WONKY AS FUCK
	# #raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $outdir/core.full.aln -p 12345 -n consensusFULL

	# #This one works for some ungodly reason. wtf and why cant i call -T for threads?
	# raxmlHPC-PTHREADS -s ${outdir}/core.full.aln -m GTRGAMMA -p 12345 -n snippy_tree

	# # bootstrap
	# raxmlHPC-PTHREADS -s core.full.aln -n snippy_tree_bootstrap -m GTRGAMMA -p 12345 -T 4 -N 100 -b 12345

	# # construct consensus MR tree
	# raxmlHPC -J MR -m GTRGAMMA -z RAxML_bootstrap.snippy_tree_bootstrap -n snippy_MR_CONS

	# # fix raxml's weird consensus tree format so topology is at least accessable
	# cat RAxML_MajorityRuleConsensusTree.snippy_MR_CONS | sed -E 's/\[[0-9]+\]//g' > fixed_snippy_MR.tre

	# sed -i -e "s/;/:0.0;/g" fixed_snippy_MR.tre

##############################################################################################################################


	#trying to clean up snippy...
	if [ $intermediate_files == "CLEAN" ]; then
	        printf "\nCleaning up intermediate output files.\n"
	        for i in $(cat ${outdir}/taxa_list.txt); do
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

	cd $outdir
	##########################################################################################
	#RAPUP RUN

	# ${rapup_path}/modules/ref_producer.py -s --ref_select $selected_ref --align_file ${gon_phy_pwd}/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas --out_file ${outdir}/${update_dir}/alignment_ref.fas

	# sed -i -e 's/.fasta//g' ${outdir}/${update_dir}/alignment_ref.fas

	# ref_name=$(head -1 ${outdir}/${update_dir}/alignment_ref.fas | sed -e 's/>//g' | sed -e 's/.ref//g')

	# printf "\nref = $ref_name\n"

	# rm "${gon_phy_pwd}/$ref_name"*

	# ls $gon_phy_pwd

	mv ${gon_phy_pwd}/trimmed_reads $outdir/trimmed_reads


	###################################################################################################33
	# IMPORTANT!
	# Output all variables (including the ones generated and used by this script) into a file
	# This is to reduce the necessity of assembling sequences when basecall is the sectio under development
	( set -o posix ; set ) > $outdir/variables.txt
	${rapup_path}/modules/var_file_fixer.py \
	--var_file ${outdir}/variables.txt \
	--out_file ${outdir}/variables.txt

	if [[ $bootstrapping == "ON" ]];
	then
		# # BOOTSTRAP RAPUP AND GONPHYLING FOR CONSENSUS TREE PRODUCTION
		# cd $outdir/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/ 
		# raxmlHPC -J MR -m GTRGAMMA -z RAxML_bootstrap.core_genome_run.out -n gon_phy_MR_CONS
		# cat RAxML_MajorityRuleConsensusTree.gon_phy_MR_CONS | sed -E 's/\[[0-9]+\]//g' > fixed_gon_phy_MR.tre
		# sed -i -e "s/;/:0.0;/g" fixed_gon_phy_MR.tre
		# sed -i -e 's/.ref//g' fixed_gon_phy_MR.tre
		# cd $outdir/rapup_run/combine_and_infer/
		# raxmlHPC -J MR -m GTRGAMMA -z RAxML_bootstrap.consensusFULL_bootstrap -n rapup_MR_CONS
		# cat RAxML_MajorityRuleConsensusTree.rapup_MR_CONS | sed -E 's/\[[0-9]+\]//g' > fixed_rapup_MR.tre
		# sed -i -e "s/;/:0.0;/g" fixed_rapup_MR.tre
		printf "\n$program_path\n"

	elif [[ $bootstrapping == "OFF" ]];
	then
		echo "BOOTSTRAPPING TURNED OFF"
	fi


elif [[ $assembly == "OFF" ]];
then
	echo "Skipping assembly and read alignment stage as per user request"
	echo "Reading variable file and beginning basecall stage"
fi




cd $outdir

##############################################################################################
# Seperating alignments into loci (if applicable) and individual taxa

if [[ $basecall == "ON" ]]; then

	source $outdir/variables.txt

	rm -r $outdir/$rapup_basecall
	rm -r $outdir/$snippy_basecall
	rm -r $outdir/$gon_phy_basecall
	rm -r $outdir/$snippy_rapup_basecall
	rm -r $outdir/gon_phy_to_rapup
	rm -r $outdir/gon_phy_to_snippy
	rm -r $outdir/rapup_to_snippy
	# rm -r $outdir/gon_phy_to_rapup/align_files
	# rm -r $outdir/gon_phy_to_snippy/align_files
	# rm -r $outdir/rapup_to_snippy/align_files
	# rm -r $outdir/gon_phy_to_rapup/assessment_output
	# rm -r $outdir/gon_phy_to_snippy/assessment_output
	# rm -r $outdir/rapup_to_snippy/assessment_output




	mkdir $outdir/$rapup_basecall
	mkdir $outdir/$snippy_basecall
	mkdir $outdir/$gon_phy_basecall
	mkdir $outdir/$snippy_rapup_basecall

	mkdir $outdir/gon_phy_to_rapup
	mkdir $outdir/gon_phy_to_snippy
	mkdir $outdir/rapup_to_snippy

	mkdir $outdir/gon_phy_to_rapup/align_files
	mkdir $outdir/gon_phy_to_snippy/align_files
	mkdir $outdir/rapup_to_snippy/align_files

	mkdir $outdir/gon_phy_to_rapup/assessment_output
	mkdir $outdir/gon_phy_to_snippy/assessment_output
	mkdir $outdir/rapup_to_snippy/assessment_output

	mkdir $outdir/gon_phy_to_snippy/blast_results
	mkdir $outdir/rapup_to_snippy/blast_results

	mkdir $outdir/$rapup_basecall/basecall_results
	mkdir $outdir/$snippy_basecall/basecall_results
	mkdir $outdir/$gon_phy_basecall/basecall_results
	mkdir $outdir/$snippy_rapup_basecall/basecall_results


	${rapup_path}/modules/locus_position_identifier.py --out_file_dir ${outdir}/${rapup_basecall}/sep_loci --position_csv_file ${outdir}/rapup_run/rapup_loci_positions.csv --concatenated_fasta ${outdir}/rapup_run/combine_and_infer/extended.aln
	# ${rapup_path}/modules/locus_position_identifier.py --out_file_dir ${outdir}/${snippy_basecall}/sep_loci --position_csv_file $loci_positions --concatenated_fasta ${outdir}/core.full.aln
	${rapup_path}/modules/locus_position_identifier.py --out_file_dir ${outdir}/${gon_phy_basecall}/sep_loci --position_csv_file $loci_positions --concatenated_fasta ${outdir}/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas

	mkdir $outdir/$snippy_basecall/sep_loci
	mkdir $outdir/$rapup_basecall/sep_loci/individual_tax
	mkdir $outdir/$snippy_basecall/sep_loci/individual_tax
	mkdir $outdir/$gon_phy_basecall/sep_loci/individual_tax

	mkdir ${outdir}/${gon_phy_basecall}/gon_phy_basecall_align
	mkdir ${outdir}/${rapup_basecall}/rapup_basecall_align
	mkdir ${outdir}/${snippy_basecall}/snippy_basecall_align

	mkdir ${outdir}/${rapup_basecall}/assessment_output
	mkdir ${outdir}/${gon_phy_basecall}/assessment_output
	mkdir ${outdir}/${snippy_basecall}/assessment_output

	
	# SEPARATE EACH LOCUS ALIGNMENT INTO INDIVIDUAL SEQUENCES (SNIPPY SEQS STAY FULL SIZE)
	for j in $(ls ${outdir}/${gon_phy_basecall}/sep_loci/*.fasta);
	do
		for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
		do
			#printf "\n$j\n$i\n"
			msa_file=$(basename $j | sed -e 's/.fasta//g')
			#printf "\n$msa_file\n"
			${rapup_path}/modules/ref_producer.py \
			-s \
			--align_file $j \
			--out_file ${outdir}/${gon_phy_basecall}/sep_loci/individual_tax/single_tax_gon_phy-${msa_file}-${i} \
			--ref_select $i
		done
	done

	for j in $(ls ${outdir}/${rapup_basecall}/sep_loci/*.fasta);
	do
	        for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
	        do
	                #printf "\n$j\n$i\n"
	                msa_file=$(basename $j | sed -e 's/.fasta//g')
			#printf "\n$msa_file\n"
	                ${rapup_path}/modules/ref_producer.py -s --align_file $j --out_file ${outdir}/${rapup_basecall}/sep_loci/individual_tax/single_tax_rapup-${msa_file}-${i} --ref_select $i
	        done
	done

	#Create single taxa locus files for snippy output
	for j in $(ls ${outdir}/core.full.aln);
	#for j in $(ls $outdir/$snippy_basecall/sep_loci/*.fasta);
	do
	        for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
	        do
	                #printf "\n$j\n$i\n"
	                msa_file=$(basename $j | sed -e 's/.fasta//g')
			#printf "\n$msa_file\n"
	                ${rapup_path}/modules/ref_producer.py \
					-s --align_file $j \
					--out_file ${outdir}/${snippy_basecall}/sep_loci/individual_tax/single_tax_snippy-${msa_file}-${i} \
					--ref_select $i
	        done
	done
	
	${program_path}/simulated_data_comparison/sim_orientation_finder.py \
	--manipulate_seqs_folder ${outdir}/${rapup_basecall}/sep_loci/individual_tax \
	--long_seqs_folder ${ref_folder} \
	--output_dir ${outdir}/${rapup_basecall}/rapup_basecall_align

	${program_path}/simulated_data_comparison/sim_orientation_finder.py \
	--long_seqs_folder ${ref_folder} \
	--manipulate_seqs_folder ${outdir}/${gon_phy_basecall}/sep_loci/individual_tax \
	--output_dir ${outdir}/${gon_phy_basecall}/gon_phy_basecall_align
	
	${program_path}/simulated_data_comparison/sim_snippy_seq_combiner.py \
	--long_seqs_folder ${ref_folder} \
	--manipulate_seqs_folder ${outdir}/${snippy_basecall}/sep_loci/individual_tax \
	--output_dir ${outdir}/${snippy_basecall}/snippy_basecall_align

	for i in $(ls ${outdir}/${rapup_basecall}/rapup_basecall_align);
	do
		mafft \
		--thread $align_threads \
		--op 5 \
		--lexp -0.5 \
		${outdir}/${rapup_basecall}/rapup_basecall_align/${i} \
		> \
		${outdir}/${rapup_basecall}/rapup_basecall_align/aligned_${i}


		${program_path}/simulated_data_comparison/sim_snippy_gapped_align_compared.py \
		-t \
		--align_1 ${outdir}/${rapup_basecall}/rapup_basecall_align/aligned_${i} \
		--output_stub rapup_assessment_${i} \
		--output_dir ${outdir}/${rapup_basecall}/assessment_output

	done

	for i in $(ls ${outdir}/${gon_phy_basecall}/gon_phy_basecall_align);
	do
		mafft \
		--thread $align_threads \
		--op 5 \
		--lexp -0.5 \
		${outdir}/${gon_phy_basecall}/gon_phy_basecall_align/${i} \
		> \
		${outdir}/${gon_phy_basecall}/gon_phy_basecall_align/aligned_${i}


		${program_path}/simulated_data_comparison/sim_snippy_gapped_align_compared.py \
		-t \
		--align_1 ${outdir}/${gon_phy_basecall}/gon_phy_basecall_align/aligned_${i} \
		--output_stub gon_phy_assessment_${i} \
		--output_dir ${outdir}/${gon_phy_basecall}/assessment_output

	done

	for i in $(ls ${outdir}/${snippy_basecall}/snippy_basecall_align);
	do
		mafft \
		--thread $align_threads \
		--op 5 \
		--lexp -0.5 \
		${outdir}/${snippy_basecall}/snippy_basecall_align/${i} \
		> \
		${outdir}/${snippy_basecall}/snippy_basecall_align/aligned_${i}


		${program_path}/simulated_data_comparison/sim_snippy_gapped_align_compared.py \
		-t \
		--align_1 ${outdir}/${snippy_basecall}/snippy_basecall_align/aligned_${i} \
		--output_stub snippy_assessment_${i} \
		--output_dir ${outdir}/${snippy_basecall}/assessment_output

	done
	
	
	
	
	
	
	
	
	
	
	
# 	ls ${outdir}/${gon_phy_basecall}/sep_loci/*.fasta | sort -R | head -${locus_num} > ${outdir}/${gon_phy_basecall}/loci_list.txt

# 	for i in $(cat ${outdir}/${gon_phy_basecall}/loci_list.txt);
# 	do
# 		cluster_name=$(basename ${i} | sed 's/.fasta//g')
# 		printf "\n$i\n"
# 		printf "\n$(basename $i)\n"
# 		${rapup_path}/modules/ref_producer.py -r --align_file $i \
# 		--out_file ${outdir}/${gon_phy_basecall}/match_finding_loci/single_tax_gon_phy-${cluster_name}
# 	done

# 	for i in $(ls $outdir/$rapup_basecall/sep_loci/*.fasta);
# 	# for i in $(ls $outdir/$rapup_basecall/sep_loci/individual_tax);
# 	do
# 		# makeblastdb -in $outdir/$rapup_basecall/sep_loci/individual_tax/$i -dbtype nucl -parse_seqids
# 		makeblastdb -in $i -dbtype nucl -parse_seqids
# 	done

# 	# BLAST THE SELECTED GON_PHYLING LOCI AGAINST THE RAPUP ALIGNMENTS
# 	# TO FIND THE MATCHING RAPUP LOCUS
# 	for i in $(ls ${outdir}/${gon_phy_basecall}/match_finding_loci/);
# 	do
# 		gon_phy_seq=$(basename $i | sed 's/.fasta//g')
# 		gon_phy_seq_path=$(realpath $i)
# 		for j in $(ls ${outdir}/${rapup_basecall}/sep_loci/*.fasta);
# 		do
# 			echo "${i}"
# 			echo "${j}"
# 			rapup_file=$(basename $j | sed 's/.fasta//g')
# 			rapup_file_path=$(realpath ${j})
# 			echo "blastn \
# 			-db ${rapup_file_path} \
# 			-query ${outdir}/${gon_phy_basecall}/match_finding_loci/$i \
# 			-out $outdir/$gon_phy_basecall/loci_finding_results/blast_output_${rapup_file}-${gon_phy_seq}.out \
# 			-max_hsps 1 \
# 			-outfmt 5"

# 			blastn \
# 			-db $j \
# 			-query ${outdir}/${gon_phy_basecall}/match_finding_loci/$i \
# 			-out ${outdir}/${gon_phy_basecall}/loci_finding_results/blast_output_${rapup_file}-${gon_phy_seq}.out \
# 			-max_hsps 1 \
# 			-outfmt 5
# 		done
# 	done

# 	#IDENTIFY EACH GON_PHY LOCUS AND ANALYZE ALL BLAST ALIGNMENTS TO FIND THE BEST ALIGNMENT
# 		${program_path}/empirical_data_comparison/emp_blast_matcher.py \
# 		--input_folder ${outdir}/${gon_phy_basecall}/loci_finding_results \
# 		--output_file ${outdir}/${gon_phy_basecall}/matched_loci_results/locus_file_matches_
	
# 	# SEPARATE EACH LOCUS ALIGNMENT INTO INDIVIDUAL SEQUENCES (SNIPPY SEQS STAY FULL SIZE)
# 	for j in $(ls ${outdir}/${gon_phy_basecall}/sep_loci/*.fasta);
# 	do
# 		for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
# 		do
# 			#printf "\n$j\n$i\n"
# 			msa_file=$(basename $j | sed -e 's/.fasta//g')
# 			#printf "\n$msa_file\n"
# 			${rapup_path}/modules/ref_producer.py \
# 			-s \
# 			--align_file $j \
# 			--out_file ${outdir}/${gon_phy_basecall}/sep_loci/individual_tax/single_tax_gon_phy-${msa_file}-${i} \
# 			--ref_select $i
# 		done
# 	done

# 	for j in $(ls ${outdir}/${rapup_basecall}/sep_loci/*.fasta);
# 	do
# 	        for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
# 	        do
# 	                #printf "\n$j\n$i\n"
# 	                msa_file=$(basename $j | sed -e 's/.fasta//g')
# 			#printf "\n$msa_file\n"
# 	                ${rapup_path}/modules/ref_producer.py -s --align_file $j --out_file ${outdir}/${rapup_basecall}/sep_loci/individual_tax/single_tax_rapup-${msa_file}-${i} --ref_select $i
# 	        done
# 	done

# 	#Create single taxa locus files for snippy output
# 	for j in $(ls ${outdir}/core.full.aln);
# 	#for j in $(ls $outdir/$snippy_basecall/sep_loci/*.fasta);
# 	do
# 	        for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
# 	        do
# 	                #printf "\n$j\n$i\n"
# 	                msa_file=$(basename $j | sed -e 's/.fasta//g')
# 			#printf "\n$msa_file\n"
# 	                ${rapup_path}/modules/ref_producer.py \
# 					-s --align_file $j \
# 					--out_file ${outdir}/${snippy_basecall}/sep_loci/individual_tax/single_tax_snippy-${msa_file}-${i} \
# 					--ref_select $i
# 	        done
# 	done

# 	snippy_seqs=$(ls $outdir/$snippy_basecall/sep_loci/individual_tax)
# 	rapup_seqs=$(ls $outdir/$rapup_basecall/sep_loci/individual_tax)
# 	gon_phyling_seqs=$(ls $outdir/$gon_phy_basecall/sep_loci/individual_tax)

# 	# NOW MAKE THE ACTUAL ALIGNMENTS AND COMPARISONS
# 	# BETWEEN RAPUP AND GON_PHYLING FIRST
# 	for j in $(cat ${outdir}/${gon_phy_basecall}/loci_list.txt);
# 	do
# 		# locus=$(basename ${j} | sed 's/.fasta//g')
# 		msa_file=$(basename ${j} | sed -e 's/.fasta//g')
# 		ln -s ${j} ${outdir}/${gon_phy_basecall}/chosen_loci/
# 		for i in $(grep ">" $j | sed -e 's/>//g' | sed -e 's/.ref//g');
# 		do
# 			#printf "\n$j\n$i\n"
# 			# msa_file=$(basename $j | sed -e 's/.fasta//g')
# 			#printf "\n$msa_file\n"
# 			${rapup_path}/modules/ref_producer.py \
# 			-s \
# 			--align_file $j \
# 			--out_file ${outdir}/${gon_phy_basecall}/chosen_loci/individual_tax/single_tax_gon_phy-${msa_file}-${i} \
# 			--ref_select $i
# 		done
# 	done

# 	for file in $(ls -1 ${outdir}/${gon_phy_basecall}/matched_loci_results/);
# 	do
# 		gon_phy_locus=$(head -1 ${outdir}/${gon_phy_basecall}/matched_loci_results/$file | \
# 		sed 's/single_tax_gon_phy-//g' | sed 's/-//g')
		
# 		rapup_locus=$(tail -1 ${outdir}/${gon_phy_basecall}/matched_loci_results/$file | \
# 		sed 's/.fasta//g' | sed 's/-//g')
		
# 		stripped_path_rapup_file=$(basename ${rapup_locus})
# 		rapup_locus=${stripped_path_rapup_file}
# 		echo $gon_phy_locus
# 		echo $rapup_locus

# 		# PROGRAM TAKES IN THE FOLDER OF SELECTED LOCI PRODUCED BY GON_PHYLING
# 		# FOLDER OF ALL RAPUP PRODUCED LOCI
# 		# THE NAME OF THE LOCUS FOR THE GON_PHYLING DATASET
# 		# THE NAME OF THE LOCUS FOR THE RAPUP DATASET
# 		# THE OUTPUT DIRECTORY FOR THE RESULTING COMBINED SEQUENCE FILES
# 		# THE OUTPUT DIRECTORY FOR THE MATCHED RAPUP SINGLE SEQUENCE FILES FOR FURTHER ALIGNMENT TO SNIPPY SEQS
# 		${program_path}/empirical_data_comparison/emp_seq_matcher.py \
# 		--folder_1 ${outdir}/${gon_phy_basecall}/chosen_loci/individual_tax \
# 		--folder_2 ${outdir}/${rapup_basecall}/sep_loci/individual_tax \
# 		--cluster_id_1 $gon_phy_locus \
# 		--cluster_id_2 $rapup_locus \
# 		--output_dir $outdir/gon_phy_to_rapup/align_files/ \
# 		--matched_seq_output_dir ${outdir}/${rapup_basecall}/chosen_loci/individual_tax

# 		for i in $(cat $outdir/gon_phy_to_rapup/align_files/taxa_name_list.txt);
# 		do
# 			mafft --op 5 --lexp -0.5 --thread $align_threads $outdir/gon_phy_to_rapup/align_files/combined_original_${gon_phy_locus}_${rapup_locus}--${i}-- > $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-original.fasta
# 			mafft --op 5 --lexp -0.5 --thread $align_threads $outdir/gon_phy_to_rapup/align_files/combined_reverse_complement_${gon_phy_locus}_${rapup_locus}--${i}-- > $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-reverse_complement.fasta
# 			mafft --op 5 --lexp -0.5 --thread $align_threads $outdir/gon_phy_to_rapup/align_files/combined_reverse_${gon_phy_locus}_${rapup_locus}--${i}-- > $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-reverse.fasta
# 			mafft --op 5 --lexp -0.5 --thread $align_threads $outdir/gon_phy_to_rapup/align_files/combined_complement_${gon_phy_locus}_${rapup_locus}--${i}-- > $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-complement.fasta

# 			${program_path}/empirical_data_comparison/emp_align_compare.py \
# 			--align_1 $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-reverse_complement.fasta \
# 			--align_2 $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-complement.fasta \
# 			--align_3 $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-reverse.fasta \
# 			--align_4 $outdir/gon_phy_to_rapup/align_files/aligned_combine-${gon_phy_locus}-${rapup_locus}--${i}-original.fasta \
# 			--output_stub $outdir/gon_phy_to_rapup/assessment_output/gon_rap_results-${gon_phy_locus}-${rapup_locus}--${i}

# 		done
# 		wait

# 	done

# 	# for i in $(ls ${outdir}/${snippy_basecall}/sep_loci/individual_tax/*);
# 	# do
# 	# 	# echo ${i}
# 	# 	makeblastdb -in $i -dbtype nucl -parse_seqids
# 	# 	snippy_genome_taxon=$(basename ${i} | sed -e 's/single_tax_snippy-core.full.aln-//g')
# 	# 	# echo ${snippy_genome_taxon}
		
		
# 	# 	for j in $(cat ${outdir}/${gon_phy_basecall}/loci_list.txt);
# 	# 	do
# 	# 		locus=$(basename ${j} | sed -e 's/.fasta//g')
			
# 	# 		for k in $(ls ${outdir}/${gon_phy_basecall}/chosen_loci/individual_tax/*${locus}*);
# 	# 		do
# 	# 			gon_phy_single_locus_single_tax=$(basename ${k} | rev | cut -d "-" -f1 | rev)
# 	# 			if [[ ${gon_phy_single_locus_single_tax} == ${snippy_genome_taxon} ]];
# 	# 			then
# 	# 				blastn \
# 	# 				-db $i \
# 	# 				-query ${k} \
# 	# 				-out $outdir/gon_phy_to_snippy/blast_results/blast_output_${locus}-${snippy_genome_taxon}.out \
# 	# 				-max_hsps 1 \
# 	# 				-outfmt 5

# 	# 			fi

# 	# 		done	
			
# 	# 	done
# 	# 	wait
# 	# done

# 	for j in $(cat ${outdir}/${gon_phy_basecall}/loci_list.txt);
# 	do
# 		locus=$(basename ${j} | sed -e 's/.fasta//g')
# 		${program_path}/empirical_data_comparison/blast_location_finder.py \
# 		--long_seqs_folder ${outdir}/snippy_basecall/sep_loci/individual_tax \
# 		--manipulate_seqs_folder ${outdir}/gon_phy_basecall/chosen_loci/individual_tax \
# 		--output_dir ${outdir}/gon_phy_to_snippy/align_files

# 	done

# 	for i in $(ls ${outdir}/gon_phy_to_snippy/align_files);
# 	do
# 		mafft \
# 		--thread $align_threads \
# 		--op 5 \
# 		--lexp -0.5 \
# 		${outdir}/gon_phy_to_snippy/align_files/${i} \
# 		> \
# 		${outdir}/gon_phy_to_snippy/align_files/aligned_${i}


# 		${program_path}/empirical_data_comparison/emp_snippy_gapped_align_compared.py \
# 		-t \
# 		--align_1 ${outdir}/gon_phy_to_snippy/align_files/aligned_${i} \
# 		--output_stub snip_gon_assessment_${i} \
# 		--output_dir ${outdir}/gon_phy_to_snippy/assessment_output

# 	done

# #####################################################################################################
# # RAPUP TO SNIPPY SECTION

# 	for j in $(ls -1 ${outdir}/${rapup_basecall}/chosen_loci/individual_tax/);
# 	do
# 		locus=$(basename ${j} | sed -e 's/single_tax-//g')
# 		taxon=$(basename ${j} | rev | cut -d "-" -f1 | rev)
# 		${program_path}/empirical_data_comparison/blast_location_finder.py \
# 		--long_seqs_folder ${outdir}/snippy_basecall/sep_loci/individual_tax \
# 		--manipulate_seqs_folder ${outdir}/rapup_basecall/chosen_loci/individual_tax \
# 		--output_dir ${outdir}/rapup_to_snippy/align_files

# 	done

# 	for i in $(ls ${outdir}/rapup_to_snippy/align_files);
# 	do
# 		mafft \
# 		--thread $align_threads \
# 		--op 5 \
# 		--lexp -0.5 \
# 		${outdir}/rapup_to_snippy/align_files/${i} \
# 		> \
# 		${outdir}/rapup_to_snippy/align_files/aligned_${i}


# 		${program_path}/empirical_data_comparison/emp_snippy_gapped_align_compared.py \
# 		-t \
# 		--align_1 ${outdir}/rapup_to_snippy/align_files/aligned_${i} \
# 		--output_stub snip_rap_assessment_${i} \
# 		--output_dir ${outdir}/rapup_to_snippy/assessment_output
# 	done

# 	#####################################################################################################
# 	# USE SINGLE SEQ FILES FROM GON_PHY AND RAPUP TO ALIGN TO SNIPPY SEQS

# # 	$program_path/empirical_data_comparison/emp_snippy_seq_matcher.py \
# # 	--manipulate_seqs_folder ${outdir}/${rapup_basecall}/chosen_loci/individual_tax \
# # 	--ref_seqs_folder $outdir/$snippy_basecall/sep_loci/individual_tax \
# # 	--align_output_dir ${outdir}/rapup_to_snippy/align_files \
# # 	--list_output_dir ${outdir}/rapup_to_snippy \
# # 	--output_align_stub snip_rap_


# # 	$program_path/empirical_data_comparison/emp_snippy_seq_matcher.py \
# # 	--manipulate_seqs_folder ${outdir}/${gon_phy_basecall}/chosen_loci/individual_tax \
# # 	--ref_seqs_folder $outdir/$snippy_basecall/sep_loci/individual_tax \
# # 	--align_output_dir ${outdir}/gon_phy_to_snippy/align_files \
# # 	--list_output_dir ${outdir}/gon_phy_to_snippy \
# # 	--output_align_stub snip_gon_

	
# # 	for tax_name in $(cat ${outdir}/gon_phy_to_snippy/taxa_list.txt);
# # 	do
# # 		for locus in $(cat ${outdir}/gon_phy_to_snippy/loci_list.txt);
# # 		do
# # 			mafft --thread $align_threads ${outdir}/gon_phy_to_snippy/align_files/snip_gon__${locus}--${tax_name}--complement.fasta > ${outdir}/gon_phy_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-complement.fasta
# # 			mafft --thread $align_threads ${outdir}/gon_phy_to_snippy/align_files/snip_gon__${locus}--${tax_name}--reverse_complement.fasta > ${outdir}/gon_phy_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-reverse_complement.fasta
# # 			mafft --thread $align_threads ${outdir}/gon_phy_to_snippy/align_files/snip_gon__${locus}--${tax_name}--reverse.fasta > ${outdir}/gon_phy_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-reverse.fasta
# # 			mafft --thread $align_threads ${outdir}/gon_phy_to_snippy/align_files/snip_gon__${locus}--${tax_name}--original.fasta > ${outdir}/gon_phy_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-original.fasta

# # 			${program_path}/empirical_data_comparison/emp_align_compare.py \
# # 			--align_1 ${outdir}/gon_phy_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-reverse_complement.fasta \
# # 			--align_2 ${outdir}/gon_phy_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-complement.fasta \
# # 			--align_3 ${outdir}/gon_phy_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-reverse.fasta \
# # 			--align_4 ${outdir}/gon_phy_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-original.fasta \
# # 			--output_stub ${outdir}/gon_phy_to_snippy/assessment_output/basecall_results-${locus}--${tax_name}-.txt
# # 		done
# # 		wait
# # 	done

# # for tax_name in $(cat ${outdir}/rapup_to_snippy/taxa_list.txt);
# # 	do
# # 		for locus in $(cat ${outdir}/rapup_to_snippy/loci_list.txt);
# # 		do
# # 			mafft --thread $align_threads ${outdir}/rapup_to_snippy/align_files/snip_rap__${locus}--${tax_name}--complement.fasta > ${outdir}/rapup_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-complement.fasta
# # 			mafft --thread $align_threads ${outdir}/rapup_to_snippy/align_files/snip_rap__${locus}--${tax_name}--reverse_complement.fasta > ${outdir}/rapup_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-reverse_complement.fasta
# # 			mafft --thread $align_threads ${outdir}/rapup_to_snippy/align_files/snip_rap__${locus}--${tax_name}--reverse.fasta > ${outdir}/rapup_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-reverse.fasta
# # 			mafft --thread $align_threads ${outdir}/rapup_to_snippy/align_files/snip_rap__${locus}--${tax_name}--original.fasta > ${outdir}/rapup_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-original.fasta

# # 			${program_path}/empirical_data_comparison/emp_align_compare.py \
# # 			--align_1 ${outdir}/rapup_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-reverse_complement.fasta \
# # 			--align_2 ${outdir}/rapup_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-complement.fasta \
# # 			--align_3 ${outdir}/rapup_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-reverse.fasta \
# # 			--align_4 ${outdir}/rapup_to_snippy/align_files/aligned_combine-${locus}--${tax_name}-original.fasta \
# # 			--output_stub ${outdir}/rapup_to_snippy/assessment_output/basecall_results-${locus}--${tax_name}-.txt
# # 		done
# # 		wait
# # 	done














# 	# # MAKE BLAST DB FOR RAPUP AND SNIPPY SEQUENCES
# 	# for i in $(ls $outdir/$snippy_basecall/sep_loci/individual_tax);
# 	# do
# 	# 	makeblastdb -in $outdir/$snippy_basecall/sep_loci/individual_tax/$i -dbtype nucl -parse_seqids
# 	# done

# 	# for i in $(ls $outdir/$rapup_basecall/sep_loci/*.fasta);
# 	# # for i in $(ls $outdir/$rapup_basecall/sep_loci/individual_tax);
# 	# do
# 	# 	# makeblastdb -in $outdir/$rapup_basecall/sep_loci/individual_tax/$i -dbtype nucl -parse_seqids
# 	# 	makeblastdb -in $outdir/$rapup_basecall/sep_loci/$i -dbtype nucl -parse_seqids
# 	# done

# 	# printf "\nBeginning actual blasting to find matching seqs regardless of names\n"
# 	# printf "\n$gon_phyling_seqs\n"
	
# 	# # MAKE LIST OF LOCI FOR ALIGNMENT AND BASECALL ASSESSMENT
# 	# ls ${gon_phyling_seqs} | sort -R | head -20 > ${outdir}/taxa_list.txt

# 	# for i in ${gon_phyling_seqs};
# 	# do
# 	# 	printf "\n$i\n"
# 	# 	blastn \
# 	# 	-db $j.fasta \
# 	# 	-query $outdir/$snippy_basecall/sep_loci/individual_tax/$i \
# 	# 	-out $outdir/$snippy_basecall/blast_results/blast_output_$cluster_grep-$file_name-$seq_grep.out \
# 	# 	-max_hsps 1 \
# 	# 	-outfmt 5
# 	# done


elif [[ $basecall == "OFF" ]]; then
printf "\nSKIPPING BASECALL SECTION OF PIPELINE PER USER SETTING\n"

fi
