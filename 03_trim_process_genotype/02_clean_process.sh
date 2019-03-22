#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N zost_clean
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=24G
#$ -t 1-17

module load intel java bwa samtools

