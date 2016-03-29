#!/bin/sh
#$ -N paq_k2_4
#$ -S /bin/sh
#$ -o p_4_out
#$ -e p_4_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-50

Rscript simulate_fix_paq.R 4
