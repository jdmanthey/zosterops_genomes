#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N var_filter
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=06:00:00
#$ -l h_vmem=8G
#$ -t 1:32

input_array=$( head -n${SGE_TASK_ID} vcf_list.txt | tail -n1 )

i=${input_array%.g.vcf}.recode.vcf


# retain only sites with variants and a minimum of 1 minor allele count
# keep only bi- or tri-allelic states

j=${i%.recode.vcf}_filtered_snps

vcftools --vcf $i --min-alleles 2 --max-alleles 3 --remove-indels --mac 1 --recode --recode-INFO-all --out $j


# process original again in different way including minimizing missing data and removing invariant sites
# for ABBA / BABA tests
# maximum number of alleles = 2 (2 alleles used in ABBA / BABA test)
# at least one sampled minor allele

j=${i%.recode.vcf}_abba

vcftools --vcf $i --min-alleles 2 --max-alleles 2 --mac 1 --max-missing 0.7 --remove-indels --recode --recode-INFO-all --out $j


# process original for demography tests 
# remove indels and need a minimum coverage of 8 per site

j=${i%.recode.vcf}_demography

vcftools --vcf $i --minDP 8 --remove-indels --recode --recode-INFO-all --out $j


