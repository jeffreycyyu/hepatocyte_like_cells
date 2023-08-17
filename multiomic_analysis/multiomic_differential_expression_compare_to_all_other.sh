#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=400G      # memory; default unit is megabytes
#SBATCH --time=0-18:00         # time (DD-HH:MM)
#SBATCH --output=multiomic_differential_expression_compare_to_all_other.out


Rscript multiomic_differential_expression_compare_to_all_other.R
