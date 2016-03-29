#!/bin/sh
#$ -N paq_boot
#$ -S /bin/sh
#$ -o p_out
#$ -e p_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-20

Rscript cluster_bootstrap_paq.R
