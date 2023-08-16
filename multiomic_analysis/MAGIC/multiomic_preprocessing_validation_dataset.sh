#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=120G      # memory; default unit is megabytes
#SBATCH --time=0-6:00         # time (DD-HH:MM)
#SBATCH --output=multiomic_preprocessing_validation_dataset.out


Rscript multiomic_preprocessing_validation_dataset.R
