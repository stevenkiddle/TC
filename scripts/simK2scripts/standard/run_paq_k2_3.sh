#!/bin/sh
#$ -N paq_k2_3
#$ -S /bin/sh
#$ -o p_3_out
#$ -e p_3_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 31-50

Rscript simulate_new_paq.R 3
