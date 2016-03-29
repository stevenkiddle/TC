#!/bin/sh
#$ -N orig
#$ -S /bin/sh
#$ -o o_1_out
#$ -e o_1_err
#$ -cwd
#$ -l h_vmem=2G
#$ -V

Rscript k2.R
