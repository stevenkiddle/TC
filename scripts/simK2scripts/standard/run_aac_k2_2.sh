#!/bin/sh
#$ -N aac_k2_2
#$ -S /bin/sh
#$ -o a_2_out
#$ -e a_2_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 31-50

Rscript simulate_new_aac.R 2
