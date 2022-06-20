#!/bin/bash

#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad
#SBATCH --job-name=HierMultinom
#SBATCH --mem-per-cpu=128gb
#SBATCH -t 4320
#SBATCH --output=logs/application_prepare_data.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

export OMP_NUM_THREADS=1

ml R/4.0

Rscript scripts/application_prepare_data.R