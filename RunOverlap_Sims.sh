#!/bin/bash

#SBATCH -o ResultsOverlap/Replicate_%a.Rout
#SBATCH --array=1-2400
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 72:00:00

module load R/3.6

R CMD BATCH --vanilla MainOverlap.R  ResultsOverlap/Replicate_${SLURM_ARRAY_TASK_ID}.Rout
