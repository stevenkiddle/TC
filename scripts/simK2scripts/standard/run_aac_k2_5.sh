#!/bin/sh
#$ -N aac_k2_5
#$ -S /bin/sh
#$ -o a_5_out
#$ -e a_5_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 31-50

Rscript simulate_new_aac.R 5
