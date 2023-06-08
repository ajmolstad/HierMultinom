#!/bin/bash

#SBATCH -o Results/Replicate_%a.Rout
#SBATCH --array=1-2400%100
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 24:00:00

module load R

R CMD BATCH --vanilla Main.R  Results/Replicate_${SLURM_ARRAY_TASK_ID}.Rout
