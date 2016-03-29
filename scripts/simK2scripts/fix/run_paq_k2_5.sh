#!/bin/sh
#$ -N paq_k2_5
#$ -S /bin/sh
#$ -o p_5_out
#$ -e p_5_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-50

Rscript simulate_fix_paq.R 5
