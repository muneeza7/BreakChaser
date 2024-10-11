#!/bin/bash

################################################################################################################
## Filename: BreakHunter.sh
## Created: June 13, 2024
## Author(s): Muneeza Maqsood, Chung Kok
##
## Purpose: 
##      - to execute BreakHunter program with conversion format implemented
##          
## Instructions:  
##      - perform all soft_clipped_reads.tsv files in the tmp_bams directory
##      - output large deletion mediated breakpoints
##
## Time estimation:
##      - 
##      - 
##         
################################################################################################################

############ parameter to set #################################
## specify output parent and tmp directory
if [ ! -d "./bams/tmp_bams/bp" ] 
then
    mkdir -pv ./bams/tmp_bams/bp
    echo "created bp folder in tmp_bams directory"
fi

if [ ! -d "./bams/tmp_bams/bp/tmp" ] 
then
    mkdir -pv ./bams/tmp_bams/bp/tmp
    echo "created tmp folder in bp directory"
fi

## pre-determined files and folders
bamsdir=./bams/tmp_bams
stage1=./script/BreakHunter_p1.R
stage2=./script/BreakHunter_p2.R
stage3=./script/PopulationSummaryBreaks.R
tmpDir=./bams/tmp_bams/bp/tmp

## extracting soft_clipped reads from bams
cd $bamsdir
echo -e "\n extracting soft clipped reads from bams... \n"

if [ $(ls | grep -c "reads.tsv$") == 0 ]; then for i in $(ls *.b37.bam); do id=${i%%.bam}; echo $id; samtools view -h $i | awk '{ if($6 ~ /S/) print $0; }' | awk '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > "${id}_soft_clipped_reads.tsv"; done; else echo -e "\n soft_clipped_reads.tsv already exist... \n"; fi

cd ../../

## BreakHunter initiating 
echo -e "\n Stage 1 initiating ... \n "
Rscript --vanilla $stage1

## Format conversion
cd $tmpDir

### making mapped.bed files:
for i in $(ls *.bed); do id=${i%%.bed}; echo "processing... $id"; bedtools intersect -a $i -b ../../../../database/RefSeq_ref_with_Intron_v1.bed -wo >X__${id}_mapped.bed ; done

### recovering intergenic breakpoints that got filtered out from mapped.bed files
for i in $(ls *_5p.bed); do id=${i%%.bed}; echo "processing... $id"; bedtools subtract -a $i -b X__${id}_mapped.bed -wa -wb > ${id}_intergenic_5p_mapped.bed; done

cd ../../../../

## stage 2
echo -e " Stage 2 initiating ... \n "
Rscript --vanilla $stage2

## stage 3
echo -e " Stage 3 initiating - Population Summary of Breaks... \n "
Rscript --vanilla $stage3

## removing tmp folder
echo -e " Files cleanning ... \n "
rm -rf $tmpDir

echo -e "\n Analysis ... completed"
