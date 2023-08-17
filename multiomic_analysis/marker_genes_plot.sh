#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=140G      # memory; default unit is megabytes
#SBATCH --time=0-3:00         # time (DD-HH:MM)
#SBATCH --output=marker_genes_plot.out


Rscript marker_genes_plot.R
