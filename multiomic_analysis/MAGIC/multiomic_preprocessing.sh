#!/bin/bash
#SBATCH --account=def-rsladek   # replace this with your own account
#SBATCH --mem-per-cpu=240G      # memory; default unit is megabytes
#SBATCH --time=0-18:00         # time (DD-HH:MM)
#SBATCH --output=multiomic_preprocessing.out


Rscript multiomic_preprocessing.R
