#! /bin/bash
# combining files for alignment and checking of snps from phycorder and tree to reads outputs

PHY_SNP_COMPARE=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

WD=$(pwd)
while getopts ":q:o:w:n:h" opt; do
  case $opt in
    q) query_align="$OPTARG"
    ;;
    o) original_fasta="$OPTARG"
    ;;
    w) working_dir="$OPTARG"
    ;;
    n) tax_num="$OPTARG"
    ;;
    h) echo  "query alignment fasta (-q), original fasta file (-o), working directory with some files and for output (-w), number of the taxon for mapping to itself in the query fasta alignment (-n)"
    exit
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

mkdir $working_dir

cd $working_dir

# pull out all sequences belonging to a particular taxon and split them into individual files
$PHY_SNP_COMPARE/alignment_splitter.py -p --align_file $query_align --taxon_num $tax_num

# cat all chunks together
cat *chunk.fas > "taxon_$tax_num-complete_phycorder_seqs.fas"

# make blastdb for the genome fasta file
makeblastdb -in $original_fasta -dbtype nucl -parse_seqids

# blast the phycorder query sequences against the original Tree to reads genome
blastn -db $original_fasta -query "taxon_$tax_num-complete_phycorder_seqs.fas" -max_target_seqs 1 -max_hsps 1 -outfmt 5 -out phycorder_to_TTR_blast_results.out

# seperate the sequences in the tree to reads genome into file
$PHY_SNP_COMPARE/alignment_splitter.py -f --align_file $original_fasta

grep -o -P "<Hsp_gaps>.{1}" phycorder_to_TTR_blast_results.out > taxon-$tax_num-gap-count.txt

grep -o "<Hsp_gaps>" phycorder_to_TTR_blast_results.out | wc -l

$PHY_SNP_COMPARE/blast_parser_fasta_matcher.py --blast_xml_file phycorder_to_TTR_blast_results.out --dir $working_dir


# TEMPFILE=/tmp/$$.tmp
# echo 0 > $TEMPFILE
# for matched_seq_file in $(ls $working_dir/matched*.txt); do
#
#   COUNTER=$[$(cat $TEMPFILE) + 1]
#
#   cat $(cat $matched_seq_file) > "unaligned_matched_seqs-$COUNTER-.fa"
#
# # TEMPFILE=/tmp/$$.tmp
# # echo 0 > $TEMPFILE
# #
# # for file in $(ls *.fas); do
# #
# #   COUNTER=$[$(cat $TEMPFILE) + 1]
# #
# #   cat $file $original_fasta > loci_master_seq_combo-$COUNTER-.fasta
# #
#   echo "$COUNTER"
#
#   echo $COUNTER > $TEMPFILE
#
# done

# unlink $TEMPFILE
