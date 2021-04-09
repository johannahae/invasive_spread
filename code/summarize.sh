#!/bin/bash
#$ -S /bin/bash
#$ -N Ksummarize.sh
#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -j y
#$ -cwd
#$ -binding linear:1
#
# Run single task in foreground
module load R/3.5.1-1
INPUT_DIR="../data/results/K/"
export INPUT_DIR
Rscript make_summary.R $INPUT_DIR 
Rscript analyse_summary.R $INPUT_DIR
Rscript check_webs.R $INPUT_DIR
#
# Script ends here
