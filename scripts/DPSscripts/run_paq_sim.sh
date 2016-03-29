#!/bin/sh
#$ -N paq_sim
#$ -S /bin/sh
#$ -o ps_out
#$ -e ps_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-50

Rscript simulate_real_tp_paq.R
