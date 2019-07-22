#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --job-name cmpfeatures
#SBATCH --output=cmp_%j.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python cmp_features.py
