#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N zost_raxml
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=01:00:00
#$ -l h_vmem=8G
#$ -t 1-19975

raxmlHPC-PTHREADS-SSE3 -T 1 -f a -x 50 -m GTRGAMMA -p 253 -N 100 \
-s /lustre/scratch/jmanthey/10_zosterops/zosterops_fasta/zost_${SGE_TASK_ID}.fasta \
-n zost_${SGE_TASK_ID}.tre -w /lustre/scratch/jmanthey/10_zosterops/zosterops_fasta/

rm /lustre/scratch/jmanthey/10_zosterops/zosterops_fasta/RAxML_bestTree.zost_${SGE_TASK_ID}.tre

rm /lustre/scratch/jmanthey/10_zosterops/zosterops_fasta/RAxML_bipartitionsBranchLabels.zost_${SGE_TASK_ID}.tre

rm /lustre/scratch/jmanthey/10_zosterops/zosterops_fasta/RAxML_bootstrap.zost_${SGE_TASK_ID}.tre

rm /lustre/scratch/jmanthey/10_zosterops/zosterops_fasta/RAxML_info.zost_${SGE_TASK_ID}.tre
