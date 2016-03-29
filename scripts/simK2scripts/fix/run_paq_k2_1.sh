#!/bin/sh
#$ -N paq_k2_1
#$ -S /bin/sh
#$ -o p_1_out
#$ -e p_1_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-50

Rscript simulate_fix_paq.R 1
