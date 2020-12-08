#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N ini_filter
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=06:00:00
#$ -l h_vmem=8G
#$ -t 1:32

input_array=$( head -n${SGE_TASK_ID} vcf_list.txt | tail -n1 )

# remove sites missing all individuals (usually low quality alignment)
# minimum depth of sequencing = 5 (liberal to begin with)
# max mean depth across individuals based on coverage distributions
# minimum quality of calls and genotype qualities of 20 

j=${input_array%.g.vcf}

vcftools --vcf $input_array --max-missing 0.05 --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 30 --recode --recode-INFO-all --out $j

