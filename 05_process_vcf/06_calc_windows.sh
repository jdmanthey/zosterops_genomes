#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N stat_windows
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1:20296

input_array=$( head -n${SGE_TASK_ID} window_list.txt | tail -n1 )

module load intel R

Rscript calculate_windows.r /lustre/scratch/jmanthey/10_zosterops/03_vcf/windows/${input_array} window_popmap.txt
