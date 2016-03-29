#!/bin/sh
#$ -N aac_k2_3
#$ -S /bin/sh
#$ -o a_3_out
#$ -e a_3_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-50

Rscript simulate_fix_aac.R 3

