# Reproduce simulation studies

**Note.** Please direct any questions to amolstad@umn.edu. 

In order to reproduce the simulation studies, one will need to first download the HierMultinom R package. 

The simulations were performed on HiperGator 3.0 at the Univeristy of Florida, which uses slurm for job management. To run the simulations with nonoverlapping coarse categories, for example, one need only execute the Run_Sims.sh script using sbatch in slurm. 

To run individual replicates locally requires only replacing the line
```
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
```
with 
``` 
uu <- a
```
where a is the replicate you seek to reproduce.  Of course, all file paths should be modified to your local system as currently, they are setup for the first author's computing account on HiperGator 3.0. 

For complete result files, which are too large to be stored on GitHub, please email the first author at the email address above. 
