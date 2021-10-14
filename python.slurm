#!/bin/bash

#SBATCH -J dskflx         # Job name
#SBATCH -o dskflx.o%j       # Name of stdout output file
#SBATCH --account=TG-AST090040
#SBATCH --partition=compute
#SBATCH -N 1                # Total # of nodes
#SBATCH --ntasks-per-node 64
#SBATCH --mem=200G
#SBATCH -t 8:00:00        # Run time (hh:mm:ss)

module load parallel

# Switch to node scratch
cd /scratch/kopec/job_$SLURM_JOB_ID/

for dir in linrot norot tctff5 tctff20 cflow fid
do

    # Extract data out of project
    echo Exracting data for $dir	
    mkdir $dir
    parallel -j64 tar -xzf {} -C $dir/ ::: $PROJECT/$dir/*tar.gz
    cd $dir
    ls .
    
    #srun python ~/enzo_storage_tools/tar_snapshots.py tar_list.txt
    #srun python ~/enzo_storage_tools/verify_tar.py
    #srun python ~/CGM-sim-analysis/save_profiles.py
    #srun python ~/CGM-sim-analysis/calc_sfr.py
    #srun python ~/CGM-sim-analysis/plot_projs.py
    #srun python ~/CGM-sim-analysis/plot_phase.py
    #srun python ~/CGM-sim-analysis/plot_pres_4panel.py
    #srun python ~/CGM-sim-analysis/save_disk_mass_growth.py
    srun  python ~/CGM-sim-analysis/save_disk_mass_flow.py
    #srun python ~/CGM-sim-analysis/plot_thermo_slabs.py
    #srun -n 32 python ~/CGM-sim-analysis/sightline_analysis_subproc_manager.py
    #srun python ~/CGM-sim-analysis/plot_gal_ev.py
    #srun python ~/CGM-sim-analysis/save_weighted_med.py
    #srun python ~/CGM-sim-analysis/save_tcool_mass_dist.py

    echo copying data products to lustre scratch
    
    # Move data products back to main scratch
    mv *txt $SCRATCH/$dir
    mv *pkl $SCRATCH/$dir
    mv *npy $SCRATCH/$dir

    # Cleanup
    echo cleanup $dir
    cd /scratch/kopec/job_$SLURM_JOB_ID/
    rm -r $dir

done
