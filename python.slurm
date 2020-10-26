#!/bin/bash

#SBATCH -J python           # Job name
#SBATCH -o python.o%j       # Name of stdout output file
#SBATCH -p development      # Queue (partition) name
#SBATCH -N 1                # Total # of nodes
#SBATCH --ntasks-per-node 68
#SBATCH -t 2:00:00        # Run time (hh:mm:ss)

#for dir in 5 20 norot linrot #fid cflow
for dir in low_mass mid_mass
do

    cd $SCRATCH/isogal_lowres_fid/$dir

    ibrun python ~/CGM_sim_analysis/save_profiles.py 
    #ibrun python ~/CGM_sim_analysis/calc_sfr.py
    #ibrun python ~/CGM_sim_analysis/plot_projs.py
    #ibrun python ~/CGM_sim_analysis/plot_phase.py
    #ibrun python ~/CGM_sim_analysis/plot_thermo_slabs.py
    #ibrun python ~/enzo_storage_tools/tar_snapshots.py tar_list.txt
    #ibrun python ~/enzo_storage_tools/verify_tar.py

done
