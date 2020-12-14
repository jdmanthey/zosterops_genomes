#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N zost_rohan
#$ -q omni
#$ -pe sm 8
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1:17

source activate rohan

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jmanthey/anaconda2/envs/rohan/lib/

/lustre/work/jmanthey/rohan/src/rohan --tstv 3.904 --size 50000 --auto 01_scaffolds.txt -t 8 -o ${SGE_TASK_ID} \
/lustre/work/jmanthey/zosterops_genome/ref.fa \
/lustre/scratch/jmanthey/10_zosterops/01_bam_files/${SGE_TASK_ID}_final.bam
