#!/bin/sh
#$ -N aac_boot
#$ -S /bin/sh
#$ -o a_out
#$ -e a_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V
#$ -t 1-20

Rscript cluster_bootstrap_aac.R
