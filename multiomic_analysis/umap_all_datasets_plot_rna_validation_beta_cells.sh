#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=140G      # memory; default unit is megabytes
#SBATCH --time=0-2:00         # time (DD-HH:MM)
#SBATCH --output=umap_all_datasets_plot_rna_validation_beta_cells.out


Rscript umap_all_datasets_plot_rna_validation_beta_cells.R
