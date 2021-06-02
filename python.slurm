#!/bin/bash

#SBATCH -J massdist         # Job name
#SBATCH -o massdist.o%j       # Name of stdout output file
#SBATCH -N 2                # Total # of nodes
#SBATCH --ntasks-per-node 24
#SBATCH -t 4:00:00        # Run time (hh:mm:ss)

for dir in  fid linrot norot tctff5 tctff20 cflow 
do

cd $SCRATCH/$dir

#ibrun python ~/enzo_storage_tools/tar_snapshots.py tar_list.txt
#ibrun python ~/enzo_storage_tools/verify_tar.py
#ibrun python ~/CGM-sim-analysis/save_profiles.py
#ibrun python ~/CGM-sim-analysis/calc_sfr.py
#ibrun python ~/CGM-sim-analysis/plot_projs.py
#ibrun python ~/CGM-sim-analysis/plot_phase.py
#ibrun -n 48 python ~/CGM-sim-analysis/save_disk_mass_growth.py
#ibrun -n 48 python ~/CGM-sim-analysis/save_disk_mass_flow.py
#ibrun -n 48 python ~/CGM-sim-analysis/plot_thermo_slabs.py
#ibrun -n 24 python ~/CGM-sim-analysis/sightline_analysis_subproc_manager.py
#ibrun -n 48 python ~/CGM-sim-analysis/plot_gal_ev.py
#ibrun -n 48 python ~/CGM-sim-analysis/save_weighted_med.py
ibrun -n 48 python ~/CGM-sim-analysis/save_tcool_mass_dist.py

done
