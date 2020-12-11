#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N zost_cut
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

# break up the depth files into single column files for each individual (locations dropped)

while read -r name1 number1; do
	number2=$((number1 + 2));
  cut zosterops_coverage.txt -f $number2 > ${name1}_depth.txt;
done < popmap.txt

