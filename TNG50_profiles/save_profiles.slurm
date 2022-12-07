#!/bin/bash

#SBATCH -J TNG50data           # Job name
#SBATCH -o TNG50data.o%j       # Name of stdout output file
#SBATCH -N 1                   # Total # of nodes
#SBATCH --ntasks-per-node 128
#SBATCH --mem=490GB
#SBATCH -t 4:00:00            # Run time (hh:mm:ss)
#SBATCH -A galaxies

module unload Python
conda activate

cd $SLURM_SUBMIT_DIR

srun -n 128 python TNG50_halo_profiles.py
