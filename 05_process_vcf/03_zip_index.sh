#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N zip_ind
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=06:00:00
#$ -l h_vmem=8G
#$ -t 1:32

input_array=$( head -n${SGE_TASK_ID} vcf_list.txt | tail -n1 )

i=${input_array%.g.vcf}.recode.vcf

bgzip -c ${i} > ${i}.gz

tabix -p vcf ${i}.gz

