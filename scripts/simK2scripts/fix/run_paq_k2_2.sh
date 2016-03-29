#!/bin/sh
#$ -N paq_k2_2
#$ -S /bin/sh
#$ -o p_2_out
#$ -e p_2_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-50

Rscript simulate_fix_paq.R 2
