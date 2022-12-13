#!/bin/bash

#SBATCH -J ext_data           # Job name
#SBATCH -o ext_data.o%j       # Name of stdout output file
#SBATCH --account=TG-AST090040
#SBATCH --partition=shared
#SBATCH -N 1                  # Total # of nodes
#SBATCH --ntasks-per-node 81
#SBATCH --mem=120GB
#SBATCH -t 24:00:00            # Run time (hh:mm:ss)

module load parallel

# Switch to node scratch
cd /scratch/kopec/job_$SLURM_JOB_ID/

for dir in cflow fid linrot norot tctff5 tctff20 
do

    # Extract data out of project
    echo Exracting data for $dir	
    mkdir $dir
    parallel -j81 tar -xzf {} -C $dir/ ::: $PROJECT/$dir/*tar.gz
    cd $dir
    # ls .
    
    #srun python ~/enzo_storage_tools/tar_snapshots.py tar_list.txt
    #srun python ~/enzo_storage_tools/verify_tar.py

    python ~/CGM-sim-analysis/save_sfh.py
    srun python ~/CGM-sim-analysis/save_sf_details.py
    srun python ~/CGM-sim-analysis/save_disk_mass_growth.py
    #srun python ~/CGM-sim-analysis/save_disk_mass_flow.py
    #srun python ~/CGM-sim-analysis/save_profiles.py
    #srun python ~/CGM-sim-analysis/save_weighted_med.py
    srun python ~/CGM-sim-analysis/save_tcool_mass_dist.py

    srun -n 40 python ~/CGM-sim-analysis/sightline_analysis_subproc_manager.py

    # If previously made SFH and masses over time, copy them over
    cp ~/data_products/$dir\_sfh.txt ./sfr.txt
    cp ~/data_products/$dir\_masses_over_time.txt ./masses_over_time.txt

    srun python ~/CGM-sim-analysis/plot_thermo_slabs.py
    #srun python ~/CGM-sim-analysis/plot_thermo_slabs_subset.py
    #srun python ~/CGM-sim-analysis/plot_phase.py
    #srun python ~/CGM-sim-analysis/plot_pres_4panel.py
    #srun python ~/CGM-sim-analysis/plot_gal_ev.py

    echo Copying data products to lustre scratch
    
    # Move data products back to main scratch
    mv *txt $SCRATCH/$dir
    mv *pkl $SCRATCH/$dir
    mv *npy $SCRATCH/$dir
    mv *png $SCRATCH/$dir
    
    # Cleanup
    echo cleanup $dir
    cd /scratch/kopec/job_$SLURM_JOB_ID/
    rm -r $dir

done
