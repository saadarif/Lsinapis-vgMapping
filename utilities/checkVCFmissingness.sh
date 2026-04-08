#!/bin/bash

#A simple scriot to check missingnes per sample in a vcf file
#Adapted From: https://darencard.net/blog/2017-01-13-missing-data-proportions-vcf/

#Provide VCF file as input
FILE=$1

#TODO: print something when there is no missing dat

paste \
<(bcftools query -f '[%SAMPLE\t]\n' ${FILE} | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' ${FILE} | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1;  count[i]++;  } END {
    for (i=1; i<=NF; i++) { 
        printf "%d\t", i;
        if (sum[i]>0)  
            print sum[i]/NR; 
        else  
            print "0%";  #WARNING: this is not working as expected
    } 
}'| sort -k1,1n | cut -f 2)
