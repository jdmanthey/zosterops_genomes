#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N zost_abba
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=01:00:00
#$ -l h_vmem=8G
#$ -t 1-1

module load intel R
Rscript abba_${SGE_TASK_ID}.r
