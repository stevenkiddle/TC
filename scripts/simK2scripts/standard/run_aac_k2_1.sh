#!/bin/sh
#$ -N aac_k2_1
#$ -S /bin/sh
#$ -o a_1_out
#$ -e a_1_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 31-50

Rscript simulate_new_aac.R 1
