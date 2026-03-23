#!/bin/bash

<<'COMMENTS'
	This is a bash version of Znolen's vgconstruct snakemake pipline:
	https://github.com/zjnolen/polyommatini-temporal-genomics/blob/main/analyses/workflow/rules/vg-construct.smk
	It first calls genotypes from modern only samples
	then uses it to create a vg-graph that can be used for read mapping
COMMENTS

##############################################
#requires: bcftools, vg
#############################################
#SETUP PARAMS
#the reference location and output location should also be set
REF=${4:-/storageToo/PROJECTS/Rebecca/WW/ref_genome/Leptidea_sinapis-GCA_905404315.1-softmasked.fa}                 
BAMFILELIST=$1 # user supplied path to all BAM files
POPLIST=$2 # a text file with individual ids in the first col and pop ids in the secon
TARGETS=$3 # REGIONS to restrict mapping too, regions excluing repeats, no sex chroms, no small scaffolds

#Depth filters
MinDP=6
MaxDP=150


#Output file conventions
OUTBCF=Modern_DTOLREF.filtered_mindp6-allsites-jointcall
#############################################

# Check if an input file is provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 BAMFILELIST POPLIST TARGETS"
    echo "See file for param details"
    exit 1
fi

#Step 1.    Calls genotypes jointly across all samples in a species dataset, excluding
#those that have been dropped from genotype call analyses. Uses bcftools
# multiallelic caller and groups individuals by sample population for the
# calling model's HWE assumption. After calling, drop low quality positions
#(QUAL < 30) and genotypes (GQ < 30, DP between 6 and 150 (~10 times average depth here). Filters to only biallelic SNPs
# with a MAF >= 0.05 and < 0.4 missing data. This ensures that only relatively
#common SNPs are used in the variation graph
bcftools mpileup -f $REF -T $TARGETS \
            -a "FORMAT/QS,FORMAT/AD,FORMAT/DP,INFO/AD" -B --min-MQ 30 \
            --min-BQ 30 -Ou --threads 10 -b $BAMFILELIST | \
            bcftools call -m -a GQ,GP -G $POPLIST -Ou --threads 10 | \
            bcftools filter -g 5 -i'QUAL >= 30' -Ou --threads 10 | \
            bcftools filter -i'FMT/GQ >= 30' -S . -Ou --threads 10 | \
            bcftools filter -i"(FMT/DP >= $MinDP) && (FMT/DP <= $MaxDP)" -S . -Ou --threads 10 | \
            bcftools view -V indels -M2 -Ou --threads 10 | \
            bcftools +fill-tags -Ou --threads 10 -- -t all | \
            bcftools filter -e'MAF<0.05 | F_MISSING>0.4' -Oz > $OUTBCF.vcf.gz
        tabix $OUTBCF.vcf.gz
	bcftools stats -s - $OUTBCF.vcf.gz > $OUTBCF.stats

#Step 2. Construct a variant graph from the reference genome and the variants
#found in the modern samples.
vg construct -r $REF -v $OUTBCF.vcf.gz -p > DTOLREF.modvars.vg

#Step 3. Indexes the variant graph.
vg index -t 10 -p -x DTOLREF.modvars.xg DTOLREF.modvars.vg
vg prune -t 10 -p -r DTOLREF.modvars.vg > DTOLREF.modvars.prune.vg
vg index -t 10 -p -g DTOLREF.modvars.gcsa DTOLREF.modvars.prune.vg

#Step 4. create snps catalog in bed format from vcf
bcftools query -f'%CHROM\t%POS0\t%END\t%REF\t%ALT\n' $OUTBCF.vcf.gz > $OUTBCF.vcf.snps.bed

