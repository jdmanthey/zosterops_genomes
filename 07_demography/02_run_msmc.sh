#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N zost_msmc
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-1

for i in $( echo Z*/ ); do
cd $i;
../msmc2_linux64bit -o $i_output -i 20 -t 6 -m 0.002166257 -p 1*2+20*1+1*2+1*3 *txt;
cd ../;
done
