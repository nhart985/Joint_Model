#!/bin/bash
#SBATCH --job-name=Simulations
#SBATCH --time=35:00
#SBATCH --mail-type=end,fail
#SBATCH --mem=500m
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2500

R CMD BATCH --no-save --no-restore Simulations.R Simulations_$SLURM_ARRAY_TASK_ID.Rout
