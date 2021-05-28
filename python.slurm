#!/bin/bash

#SBATCH -J massflx          # Job name
#SBATCH -o massflx.o%j       # Name of stdout output file
#SBATCH -N 2                # Total # of nodes
#SBATCH --ntasks-per-node 24
#SBATCH -t 3:00:00        # Run time (hh:mm:ss)

for dir in tctff20 cflow fid linrot norot tctff5 
do

cd $SCRATCH/$dir

#ibrun python ~/CGM_sim_analysis/save_profiles.py
#ibrun python ~/CGM_sim_analysis/calc_sfr.py
#ibrun python ~/CGM_sim_analysis/plot_projs.py
#ibrun python ~/CGM_sim_analysis/plot_phase.py
#ibrun -n 48 python ~/CGM_sim_analysis/save_quantities_over_time.py
ibrun -n 48 python ~/CGM_sim_analysis/save_disk_mass_flow.py
#ibrun -n 48 python ~/CGM_sim_analysis/plot_thermo_slabs.py
#ibrun python ~/enzo_storage_tools/tar_snapshots.py tar_list.txt
#ibrun python ~/enzo_storage_tools/verify_tar.py
#ibrun -n 24 python ~/CGM_sim_analysis/sightline_analysis_subproc_manager.py
#ibrun -n 48 python ~/CGM_sim_analysis/plot_gal_ev.py
#ibrun -n 48 python ~/CGM_sim_analysis/save_weighted_med.py

done
