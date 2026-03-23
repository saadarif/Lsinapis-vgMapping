#!/bin/bash

<<'COMMENTS'
	This is a bash version of Znolen's ref_filt.smk snakemake pipline:
	https://github.com/zjnolen/PopGLen/blob/master/workflow/rules/0.2_ref_filt.smk
	It takes in a fasta, fasta index, repeat bed file for thta fasta, any (sex) chromosomes
	to exclude and any contigs less than a particular size 
	then uses it to creates a mappability mask using genamp
	It provides a bed file to target sites that exclude:
	repeats
	mappable regions below a threshold
	any (sex) chromosomes or others
	any contig less than a specified size
	Filters for depth and missingness are applied via bcf tools during calling genotypes
COMMENTS

##############################################
#requires: bedtools, genmap, bedops
#############################################
#SETUP PARAMS
#Tools paths
bedtools=/opt/bedtools

#the reference location and output location should also be set
REF=/storageToo/PROJECTS/Rebecca/WW/ref_genome/Leptidea_sinapis-GCA_905404315.1-softmasked.fa                 
REFFAI=/storageToo/PROJECTS/Rebecca/WW/ref_genome/Leptidea_sinapis-GCA_905404315.1-softmasked.fa.fai

#Repeat library as bed file
REPBED=/storageToo/PROJECTS/Rebecca/WW/ref_genome/Leptidea_sinapis-GCA_905404315.1-softmasked.repeats.sorted.bed.sorted

#Collapsed regions in the DTOL assembly, provided by JESPER BOMAN
# generated using control-FREEC
COLL=/storageToo/PROJECTS/Rebecca/WW/ref_genome/Boman_FD_CAT_FD_SWE_10_at_least_two_overlaps_matchChrNames_sorted.bed 

#Exclusion filters
SEXCHR=(Z 2 3) #exclude sex chromomes in DTOL ref
MINSIZE=1000000 #exclude an contigs less than 1MB

#GENMAP filters
K=30
E=2
MAP_THRESH=1

#Output file conventions
OUTPREFIX=Leptidea_sinapis-GCA_905404315.1-softmasked
##################################################
#All files should be SORTED using sort-bed (BEDOPS)
#Step 1. Create a bed file for the entire genome
awk -v OFS='\t' '{{print $1, "0", $2}}' $REFFAI > $OUTPREFIX.bed
sort-bed $OUTPREFIX.bed > $OUTPREFIX.sorted.bed
#summarize bed
len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' $OUTPREFIX.sorted.bed)
echo $len | awk '{{print "Total genome\t"$1"\t"$1/$1*100}}' > $OUTPREFIX.bed.sum

#Step 2. Make bed files of sex chromosomes and other excluded regrions
#These will be ssubtracted to make the final targets bed file later
# Join array elements with a pipe (|) for regex: "Z|2|3"
REGEX_PATTERN=$(IFS='|'; echo "${SEXCHR[*]}")
# Filter for lines where column 1 ($1) matches the pattern exactly
awk -v pat="^($REGEX_PATTERN)$" '$1 ~ pat' $OUTPREFIX.sorted.bed > ${OUTPREFIX}_sex_chromosomes.bed
# summarize bed
len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' ${OUTPREFIX}_sex_chromosomes.bed)
echo $len | awk '{{print "Total genome\t"$1"\t"$1/$1*100}}' > ${OUTPREFIX}_sex_chromosomes.bed.sum

#Create a bed file of all scaffolds under a specified size
awk  -v min="$MINSIZE" '$3 < min' $OUTPREFIX.sorted.bed > ${OUTPREFIX}_scaffLess1MBbp.bed
len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' ${OUTPREFIX}_scaffLess1MBbp.bed)
echo $len | awk '{{print "Total genome\t"$1"\t"$1/$1*100}}' > ${OUTPREFIX}_scaffLess1MBbp.bed.sum


#Step 3. Index reference for Genmap
genmap index -F $REF -I ./genmapIndex
#Estimate mappability of each site in the genome
genmap map -K $K -E $E -I ./genmapIndex -O ${OUTPREFIX}_k${K}_e${E} -bg 
#Generate sliding windows across the genome to average mappability over.
$bedtools makewindows -b $OUTPREFIX.bed -w $K -s 1 > ${OUTPREFIX}_windows_k${K}.bed
#Generate a bed file containing the pileup mappability of a site, i.e. the
#mean mappability of all possible kmers mapping to it.
awk '{print $1"\t"$2"\t"$3"\t"$2"-"$3"\t"$4}' ${OUTPREFIX}_k${K}_e${E}.bedgraph > ${OUTPREFIX}_k${K}_e${E}.bedgraph.tmp
#sort files before bedmap
sort-bed ${OUTPREFIX}_k${K}_e${E}.bedgraph.tmp > ${OUTPREFIX}_k${K}_e${E}.bedgraph.tmp.sorted
sort-bed ${OUTPREFIX}_windows_k${K}.bed > ${OUTPREFIX}_windows_k${K}.bed.sorted

bedmap --echo --wmean  ${OUTPREFIX}_windows_k${K}.bed.sorted ${OUTPREFIX}_k${K}_e${E}.bedgraph.tmp.sorted | tr "|" "\t" | \
        awk '{print $1"\t"$3-1"\t"$3"\t"$4}' > ${OUTPREFIX}_pileup_mappability_k${K}_e${E}.bed
#cleanup 
# Cleanup intermediate files
rm "${OUTPREFIX}_k${K}_e${E}.bedgraph.tmp"*
rm "${OUTPREFIX}_windows_k${K}.bed.sorted"
# Create a bed containing all sites with a mappability below a set threshold
# These will be subtracted from the final target files
awk -v map_thresh="$MAP_THRESH" '$4 < map_thresh' ${OUTPREFIX}_pileup_mappability_k${K}_e${E}.bed \
	 > $OUTPREFIX_k${K}_e${E}_thresh${MAP_THRESH}.bed.tmp
$bedtools merge -i $OUTPREFIX_k${K}_e${E}_thresh${MAP_THRESH}.bed.tmp > $OUTPREFIX_k${K}_e${E}_thresh${MAP_THRESH}.bed
rm $OUTPREFIX_k${K}_e${E}_thresh${MAP_THRESH}.bed.tmp

#summarize bed
len=$(awk 'BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}' $OUTPREFIX_k${K}_e${E}_thresh${MAP_THRESH}.bed)
echo $len $(awk -F "\t" '{print $2}' $OUTPREFIX.bed.sum) | awk -v k="$K" -v e="$E" -v th="$MAP_THRESH" \
    '{print "Pileup mappability K"k"-E"e" <"th"\t"$2-$1"\t"($2-$1)/$2*100}' \
    > "${OUTPREFIX}_k${K}_e${E}_thresh${MAP_THRESH}.bed.sum"


#Step 4. Remove all bed generated above to provide a bed file for target sites to include during genotype calling

#remove sexchromosomes
$bedtools subtract -a ${OUTPREFIX}.sorted.bed -b ${OUTPREFIX}_sex_chromosomes.bed > ${OUTPREFIX}_autosomes.bed
#remove scaffolds less than minimum size
$bedtools subtract -a ${OUTPREFIX}_autosomes.bed -b ${OUTPREFIX}_scaffLess1MBbp.bed > ${OUTPREFIX}_autosomes_minScaf1MB.bed
#remove any repeat regions
$bedtools subtract -a  ${OUTPREFIX}_autosomes_minScaf1MB.bed -b $REPBED > ${OUTPREFIX}_autosomes_minScaf1MB_repma.bed
#remove regions below mappability threshold
$bedtools subtract -a ${OUTPREFIX}_autosomes_minScaf1MB_repma.bed -b $OUTPREFIX_k${K}_e${E}_thresh${MAP_THRESH}.bed > ${OUTPREFIX}_autosomes_minScaf1MB_repma_snpabble.bed
#remove any collapsed regions from the assembly
$bedtools subtract -a ${OUTPREFIX}_autosomes_minScaf1MB_repma_snpabble.bed -b $COLL > ${OUTPREFIX}_autosomes_minScaf1MB_repma_snpabble_collma.bed

# summarize bed
len=$(awk 'BEGIN{{SUM=0}}{{SUM+=$3-$2}}END{{print SUM}}' ${OUTPREFIX}_autosomes_minScaf1MB_repma_snpabble_collma.bed)
echo $len  $(awk -F "\t" '{print $2}' $OUTPREFIX.bed.sum)| awk '{{print "Total genome\t"$1"\t"$1/$2*100}}' > ${OUTPREFIX}_autosomes_minScaf1MB_repma_snpabble_collma.bed.sum
