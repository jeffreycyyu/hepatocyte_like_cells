#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=120G      # memory; default unit is megabytes
#SBATCH --time=0-0:30         # time (DD-HH:MM)
#SBATCH --output=proportions_tissue_cell_type_markers.out


Rscript proportions_tissue_cell_type_markers.R
