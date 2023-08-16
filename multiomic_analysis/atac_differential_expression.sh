#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=500G      # memory; default unit is megabytes
#SBATCH --time=0-26:00         # time (DD-HH:MM)
#SBATCH --output=atac_differential_expression.out


Rscript atac_differential_expression.R
