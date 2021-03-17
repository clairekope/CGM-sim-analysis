#!/usr/bin/env python
# coding: utf-8
import yt
yt.enable_parallelism()
import numpy as np
import matplotlib.pyplot as plt
import pickle

n_theta = 13
phi_step = np.deg2rad(10)
files = ['entropy','cooling_time','density','temperature','metallicity','pressure']
units = ['keV*cm**2','Gyr','g*cm**-3','K','Zsun','dyne*cm**-2']


# currently limit theta to upper hemisphere
theta = np.linspace(0, np.pi, n_theta, endpoint=True)

phi = np.arange(0, 2*np.pi, phi_step)
n_phi = phi.size
print(n_phi)

theta_coord, phi_coord = np.meshgrid(theta,phi) # each row is const phi

results = {}

datasets = yt.load("DD????/DD????")
for sto, ds in datasets.piter(dynamic=True, storage=results):

    length = ds.quan(206, 'kpc')
    s = 0.5 * np.ones((3, n_phi, n_theta)) # dim 1 is const phi; dim 0 is x,y,z
    starts = ds.arr(s, 'code_length')
    uvec = np.array([np.sin(theta_coord)*np.cos(phi_coord),
                     np.sin(theta_coord)*np.sin(phi_coord),
                     np.cos(theta_coord)])
    ends = starts + length*uvec

    # each list within rays will have constant theta
    rays = []

    # the phi angles for i_theta==0 and i_theta==n_theta-1 are all redundant
    # as these are the poles
    rays.append(ds.ray(starts[:, 0, 0], ends[:, 0, 0]))
    for i_theta in range(1,n_theta-1):
        rays.append([])
        for i_phi in range(n_phi):
            rays[i_theta].append(ds.ray(starts[:, i_phi, i_theta], ends[:, i_phi, i_theta]))
    rays.append(ds.ray(starts[:, -1, -1], ends[:, -1, -1]))

    r_edges = np.logspace(np.log10(2e-1), np.log10(206), 21)
    r_centers = r_edges[:-1] + np.diff(r_edges)/2

    # Is this the most efficient outer loop? heck no
    # It should really be the innermost loop
    # But it's easy to program
    for quantity_name, unit in zip(fields, units):
        
        quantity_arrays[quantity_name] = {}

        # Bin quantity of all rays with matching theta
        for i_theta in range(0,n_theta):

            quantities = []

            # Theta = 0 or 2pi
            if i_theta==0 or i_theta==n_theta-1:
                ray = rays[0]
                r_binner = np.digitize(ray['radius'].to('kpc'), r_edges)
                for i in range(1, r_edges.size):
                    this_bin = r_binner==i
                    quantities.append(list(ray[quantity_name][this_bin].to(unit).value))

            else:
                # first ray at this theta 
                ray = rays[i_theta][0]    
                r_binner = np.digitize(ray['radius'].to('kpc'), r_edges)
                for i in range(1, r_edges.size):
                    this_bin = r_binner==i
                    quantities.append(list(ray[quantity_name][this_bin].to(unit).value))

                # subsequent rays at this theta
                for ray in rays[i_theta][1:]:
                    r_binner = np.digitize(ray['radius'].to('kpc'), r_edges)
                    for i in range(1, r_edges.size):
                        this_bin = r_binner==i
                        quantities[i-1].extend(ray[quantity_name][this_bin].to(unit).value)

            # Process/compress bins into med, min, and max
            quantity_min = np.zeros(len(quantities))
            quantity_low = np.zeros(len(quantities))
            quantity_med = np.zeros(len(quantities))
            quantity_upp = np.zeros(len(quantities))
            quantity_max = np.zeros(len(quantities))

            for i in range(len(quantities)):
                quantity_min[i] = np.min(quantities[i])
                quantity_low[i] = np.percentile(quantities[i], 25)
                quantity_med[i] = np.median(quantities[i])
                quantity_upp[i] = np.percentile(quantities[i], 75)
                quantity_max[i] = np.max(quantities[i])

            quantity_arrays[quantity_name][i_theta] = {'min':quantity_min,
                                                       'p25':quantity_low,
                                                       'med':quantity_med,
                                                       'p75':quantity_upp,
                                                       'max':quantity_max}
            
    sto.result = quantity_arrays
    sto.result_id = ds.basename

if yt.is_root():
    pickle.dump(results, "hedgehog_data.pkl", protocol=3)
