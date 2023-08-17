#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=140G      # memory; default unit is megabytes
#SBATCH --time=0-1:00         # time (DD-HH:MM)
#SBATCH --output=rna_teng2018cellreports_gene_per_cell_type_histogram_display_non_imputed_proportion_zeros.out


Rscript rna_teng2018cellreports_gene_per_cell_type_histogram_display_non_imputed_proportion_zeros.R
