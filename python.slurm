#!/bin/bash

#SBATCH -J tcool_all         # Job name
#SBATCH -o tcool_all.o%j       # Name of stdout output file
#SBATCH --account=TG-AST090040
#SBATCH --partition=shared
#SBATCH -N 1                # Total # of nodes
#SBATCH --ntasks-per-node 64
#SBATCH --mem=128G
#SBATCH -t 2:00:00        # Run time (hh:mm:ss)

module purge
module r

for dir in linrot norot tctff5 tctff20 cflow fid
do

cd $SCRATCH/$dir

#srun -n 64 python ~/enzo_storage_tools/tar_snapshots.py tar_list.txt
#srun -n 64 python ~/enzo_storage_tools/verify_tar.py
#srun -n 64 python ~/CGM-sim-analysis/save_profiles.py
#srun -n 64 python ~/CGM-sim-analysis/calc_sfr.py
#srun -n 64 python ~/CGM-sim-analysis/plot_projs.py
#srun -n 64 python ~/CGM-sim-analysis/plot_phase.py
#srun -n 64 python ~/CGM-sim-analysis/plot_pres_4panel.py
#srun -n 64 python ~/CGM-sim-analysis/save_disk_mass_growth.py
#srun -n 64 python ~/CGM-sim-analysis/save_disk_mass_flow.py
#srun -n 64 python ~/CGM-sim-analysis/plot_thermo_slabs.py
#srun -n 32 python ~/CGM-sim-analysis/sightline_analysis_subproc_manager.py
#srun -n 64 python ~/CGM-sim-analysis/plot_gal_ev.py
#srun -n 64 python ~/CGM-sim-analysis/save_weighted_med.py
srun -n 64 python ~/CGM-sim-analysis/save_tcool_mass_dist.py

done
