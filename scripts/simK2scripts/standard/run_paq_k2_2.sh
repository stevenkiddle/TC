#!/bin/sh
#$ -N paq_k2_2
#$ -S /bin/sh
#$ -o p_2_out
#$ -e p_2_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 31-50

Rscript simulate_new_paq.R 2
