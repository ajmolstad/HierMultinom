#!/bin/bash

#SBATCH -o results_new/M1_%a.Rout
#SBATCH --array=1-450
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad
#SBATCH --mem-per-cpu=12gb
#SBATCH -t 144:00:00

export OMP_NUM_THREADS=1

module load R

R CMD BATCH --vanilla application1.R  results/M1_${SLURM_ARRAY_TASK_ID}.Rout
