#!/bin/bash
#SBATCH --job-name=cg
#SBATCH --output=cg_%j.out
#SBATCH --error=cg_%j.err
#SBATCH --ntasks=4
#SBATCH --time=00-00:10:00
#SBATCH --constraint=highmem

module load impi

srun ./cg -r benchmark/4_cg_random.in