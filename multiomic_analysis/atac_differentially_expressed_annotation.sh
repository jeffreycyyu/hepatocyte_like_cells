#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=80G      # memory; default unit is megabytes
#SBATCH --time=0-1:00         # time (DD-HH:MM)
#SBATCH --output=atac_differentially_expressed_annotation.out


Rscript atac_differentially_expressed_annotation.R
