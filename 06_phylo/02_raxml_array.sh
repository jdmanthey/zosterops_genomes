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
#$ -l h_rt=20:00:00
#$ -l h_vmem=8G
#$ -t 1-20259

input_file=$( head -n${SGE_TASK_ID} helper_fasta.txt | tail -n1 )

raxmlHPC-PTHREADS-SSE3 -T 1 -f a -x 50 -m GTRGAMMA -p 253 -N 100 \
-s /lustre/scratch/jmanthey/10_zosterops/03_vcf/fasta/${input_file} \
-n ${input_file%.fasta}.tre \
-w /lustre/scratch/jmanthey/10_zosterops/03_vcf/trees/

rm /lustre/scratch/jmanthey/10_zosterops/03_vcf/trees/RAxML_bestTree.${input_file%.fasta}.tre

rm /lustre/scratch/jmanthey/10_zosterops/03_vcf/trees/RAxML_bipartitionsBranchLabels.${input_file%.fasta}.tre

rm /lustre/scratch/jmanthey/10_zosterops/03_vcf/trees/RAxML_bootstrap.${input_file%.fasta}.tre

rm /lustre/scratch/jmanthey/10_zosterops/03_vcf/trees/RAxML_info.${input_file%.fasta}.tre

