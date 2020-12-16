#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N var_simp
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=06:00:00
#$ -l h_vmem=8G
#$ -t 1:32

input_array=$( head -n${SGE_TASK_ID} vcf_list.txt | tail -n1 )

i=${input_array%.g.vcf}_demography.recode.vcf

j=${i%_demography.recode.vcf}_demography.simple.vcf

# use bcftools to simplify the vcf files to reduce file size, complexity, and make them easier to work with

bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' $i > $j
