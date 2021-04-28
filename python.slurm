#!/bin/bash

#SBATCH -J thermo           # Job name
#SBATCH -o thermo.o%j       # Name of stdout output file
#SBATCH -N 2                # Total # of nodes
#SBATCH --ntasks-per-node 24
#SBATCH -t 12:00:00        # Run time (hh:mm:ss)

for dir in fid linrot norot tctff5 tctff20 cflow
do

cd $SCRATCH/$dir

#ibrun python ~/CGM_sim_analysis/save_profiles.py
#ibrun python ~/CGM_sim_analysis/calc_sfr.py
#ibrun python ~/CGM_sim_analysis/plot_projs.py
#ibrun python ~/CGM_sim_analysis/plot_phase.py
#ibrun -n 48 python ~/CGM_sim_analysis/save_quantities_over_time.py
#ibrun -n 48 python ~/CGM_sim_analysis/plot_thermo_slabs.py
#ibrun python ~/enzo_storage_tools/tar_snapshots.py tar_list.txt
#ibrun python ~/enzo_storage_tools/verify_tar.py
#ibrun -n 24 python ~/CGM_sim_analysis/sightline_analysis_subproc_manager.py
#ibrun -n 48 python ~/CGM_sim_analysis/plot_gal_ev.py
ibrun -n 48 python ~/CGM_sim_analysis/save_weighted_med.py

done
