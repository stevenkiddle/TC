#!/bin/sh
#$ -N aac_k2_4
#$ -S /bin/sh
#$ -o a_4_out
#$ -e a_4_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-50

Rscript simulate_fix_aac.R 4
