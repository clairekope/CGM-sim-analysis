# Contents

A collection of scripts for plotting & extracting data from my CGM simulations.

| File prefix | Purpose |
|-------------|---------|
| calc        | Helper functions that can be imported by "save" files |
| gen         | Generate files for simulation manangement (no data extraction) |
| plot        | Save plots! |
| save        | Extract data from simulation outputs and save it, either as ASCII tables or Python 3 pickled dictionaries of dictionaries |
| sightline_analysis | See below |

## Misc Files

* annotate_phase_plot - semi-automatedly add lines of constant pressure & entropy to a density-temperature phase plot from yt.
* test_HSE.py - explore deviations from hydrostatic equilibrium that may be present in initial conditions
* manual_profile_snippet.txt - an exerpt of code showing how to make yt-esque profiles with SciPy (more flexible, too).

# Running on Comet (SLURM)

Use `sbatch python.slurm` to submit. Un/comment the lines according to which script you want to run.
It's low tech, man. If you're someone else using this code, consider this file a template.

# Sightline Analysis Scripts

There are a pair of scripts with this prefix. The file `sightline_analysis.py` is executed by subprocesses
launched and managed by `sightline_analysis_subproc_manager.py`. This is because the data extraction in 
`sightline_analysis.py` can 1) fail due to bad handling of the simulation grids, requiring a reload of the dataset
and 2) segfault. When running `sightline_analysis_subproc_manager.py` in parallel, ask for twice as many processes
as you would like to have handling the data. This is because half of the set will be "minding" the other half.

Use `collate_pkls.py` to combine the resulting individual pickle files into one.
