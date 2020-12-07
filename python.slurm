#!/bin/bash

#SBATCH -J python           # Job name
#SBATCH -o python.o%j       # Name of stdout output file
#SBATCH -N 1                # Total # of nodes
#SBATCH --ntasks-per-node 28
#SBATCH -t 4:00:00        # Run time (hh:mm:ss)

#for dir in lowmass midmass
#do

cd $SCRATCH/lowres-midmass

#srun python ~/CGM_sim_analysis/save_profiles.py
srun python ~/CGM_sim_analysis/calc_sfr.py
#srun python ~/CGM_sim_analysis/plot_projs.py
#srun python ~/CGM_sim_analysis/plot_phase.py
srun python ~/CGM_sim_analysis/plot_thermo_slabs.py
#srun python ~/enzo_storage_tools/tar_snapshots.py tar_list.txt
#srun python ~/enzo_storage_tools/verify_tar.py

#done
