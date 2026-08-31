#!/bin/bash
#SBATCH --job-name=Combine
#SBATCH --time=10:00
#SBATCH --mail-type=fail
#SBATCH --mem=250m
#SBATCH --cpus-per-task=1

R CMD BATCH --no-save --no-restore Combine.R Combine.Rout
