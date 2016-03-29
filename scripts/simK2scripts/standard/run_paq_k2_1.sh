#!/bin/sh
#$ -N paq_k2_1
#$ -S /bin/sh
#$ -o p_1_out
#$ -e p_1_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 31-50

Rscript simulate_new_paq.R 1
