##############################################################
# Save whole-galaxy mass-weighted medians of listed quantities
##############################################################

import yt
yt.enable_parallelism()
import numpy as np
import wquantiles as wq
import pickle

radial_bins = 51

quantities = ['cooling_time','entropy','pressure']
units = ['Gyr', 'keV*cm**2', 'dyn*cm**-2']
profiles = {}

datasets = yt.load('DD????/DD????')
for sto, ds in datasets.piter(storage=profiles, dynamic=False):
    
    sph = ds.sphere('c', (206,'kpc'))
    r_edges = np.linspace(2e-1, 206, radial_bins) * yt.units.kpc    
    r_binner = np.digitize(sph['radius'], r_edges.to('cm'))
    
    quantity_arrays = {}
    
    for quantity, unit in zip(quantities, units):

        prof_med = np.ones(r_edges.size-1) * np.nan
        prof_16 = np.ones(r_edges.size-1) * np.nan
        prof_84 = np.ones(r_edges.size-1) * np.nan

        for i in range(1, r_edges.size):
            this_bin = r_binner==i
            try:
                prof_med[i-1] = wq.median(sph[quantity][this_bin].to(unit),
                                        sph['cell_mass'][this_bin])
                prof_16[i-1] = wq.quantile(sph[quantity][this_bin].to(unit),
                                         sph['cell_mass'][this_bin],
                                         0.16)
                prof_84[i-1] = wq.quantile(sph[quantity][this_bin].to(unit),
                                         sph['cell_mass'][this_bin],
                                         0.84)
            except ValueError:
                continue

        quantity_arrays[quantity] = {'med': prof_med,
                                     'p16': prof_16,
                                     'p84': prof_84}
    
    sto.result = quantity_arrays
    sto.result_id = ds.basename
    
with open("weighted_med_profs.pkl", "wb") as f:
    pickle.dump(profiles, f, protocol=3)
