#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=400G      # memory; default unit is megabytes
#SBATCH --time=0-5:00         # time (DD-HH:MM)
#SBATCH --output=marker_genes_plot_gene_set.out


Rscript marker_genes_plot_gene_set.R
