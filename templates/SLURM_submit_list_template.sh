#!/bin/bash

#SBATCH --account=ivsusers
#SBATCH --job-name=JOBNAME
#SBATCH --time WALLTIME
#SBATCH --output=/dev/null
#SBATCH --error=./SLURM.err
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=CPU
#SBATCH --mem=MEMORY
#SBATCH --exclude=pleiad14,pleiad15

list=LIST

line=$( sed "${SLURM_ARRAY_TASK_ID}q;d" $list )

srun bash $line
