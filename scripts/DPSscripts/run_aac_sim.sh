#!/bin/sh
#$ -N aac_sim
#$ -S /bin/sh
#$ -o aac_out
#$ -e aac_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-50

Rscript simulate_real_tp_aac.R
