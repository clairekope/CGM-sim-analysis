#!/usr/bin/env python
# coding: utf-8

import yt
yt.enable_parallelism()
import pickle
import numpy as np

bins = 51

fields = [('gas','entropy'),
        ('gas','cooling_time'),
        ('gas','density'),
        ('gas','temperature'),
        ('gas','metallicity'),
        ('gas','pressure'),
        ('gas','cell_mass'),
        ('gas','radial_velocity')]

units = ['keV*cm**2','Gyr','g*cm**-3',
         'K','Zsun','dyne*cm**-2',
         'Msun','km/s']

theta = [0,]
theta.extend(np.arange(15,180,30))
theta.append(180)

edges = np.linspace(2e-1, 206, bins)
centers = edges[:-1] + np.diff(edges)/2

datasets = yt.load("DD????/DD????")
storage = {}

for my_storage, ds in datasets.piter(storage=storage):

    quantity_arrays = {}
    for quantity_name in fields:
        quantity_arrays[quantity_name[1]] = {}

    sph = ds.sphere('c',(206,'kpc'))

    for t in range(len(theta)-1):

        select = f"(obj[('index','spherical_theta')] > {np.deg2rad(theta[t])}) &"\
                 f"(obj[('index','spherical_theta')] < {np.deg2rad(theta[t+1])})"

        wedge = sph.cut_region([select])
        wedge.set_field_parameter('center',ds.arr([0.5,0.5,0.5], 'code_length'))
        wedge.get_data(('index','radius'))

        binner = np.digitize(wedge[('index','radius')].to('kpc'), edges)

        for quantity_name, unit in zip(fields, units):

            quantities = []

            for i in range(1, edges.size):
                this_bin = binner==i
                quantities.append(list(wedge[quantity_name][this_bin].to(unit).value))

            quantity_min = np.zeros(len(quantities))
            quantity_low = np.zeros(len(quantities))
            quantity_med = np.zeros(len(quantities))
            quantity_upp = np.zeros(len(quantities))
            quantity_max = np.zeros(len(quantities))

            for i in range(len(quantities)):
                quantity_min[i] = np.min(quantities[i])
                quantity_low[i] = np.percentile(quantities[i], 16)
                quantity_med[i] = np.median(quantities[i])
                quantity_upp[i] = np.percentile(quantities[i], 84)
                quantity_max[i] = np.max(quantities[i])

            quantity_arrays[quantity_name[1]][t] = {'min':quantity_min,
                                                 'p16':quantity_low,
                                                 'med':quantity_med,
                                                 'p84':quantity_upp,
                                                 'max':quantity_max}

    my_storage.result = quantity_arrays
    my_storage.result_id = ds.basename

    ds.close()
    
if yt.is_root():
    with open("wedges.pkl","wb") as f:
        pickle.dump(storage, f, protocol=3)
