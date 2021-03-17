import yt
yt.enable_parallelism()

import yt.units as u
import numpy as np
import os
from calc_enclosed_mass import *

quants = ['density','temperature','entropy','pressure',
          'metallicity', 'cooling_time']

datasets = yt.load("DD????/DD????")

if os.path.isfile('last_profile'):
    with open('last_profile','r') as f:
        last_profile = int(f.read())

else:
    last_profile = -1

for ds in datasets.piter(dynamic=True, ):
    if int(ds.basename[-4:]) > last_profile:
        center = ds.quan(0.5, 'code_length')
        le = center-ds.quan(0.6, 'Mpc')
        re = center+ds.quan(0.6, 'Mpc')
        offset = ds.quan(5, 'kpc')

        sph = ds.sphere('c',(500,'kpc'))        
        dsk_slb = ds.region([center, center, center],
                            [le, le, center-offset],
                            [re, re, center+offset])

        cgm = sph-dsk_slb
        cgm.set_field_parameter('center', ds.arr([0.5,0.5,0.5],'code_length'))

        prof = yt.create_profile(cgm, 'radius', quants, weight_field='cell_mass',
                                 units = {'radius':'kpc', 'cooling_time':'Gyr',
                                          'entropy':'keV*cm**2'})

        arr = np.empty((prof.x.size, len(quants)+2)) # add radius and tcool/tff
        arr[:,0] = prof.x
        for i, field in enumerate(quants):
            arr[:,i+1] = prof[field]

        # total NFW, gas, and particle mass. Average Miyamoto & Nagai acceleration
        Mtot_r = NFW_mass_enclosed(prof.x) + \
                 cell_and_particle_mass_enclosed(cgm, prof.x_bins.v)
    
        g_r = u.G * Mtot_r / np.power(prof.x,2)
        g_r = g_r.to('cm/s**2') + MN_accel(prof.x)
    
        arr[:,-1] = prof['cooling_time'] / np.sqrt(2*prof.x/g_r)
    
        header = 'radius'
        for field in quants:
            header += ' ' + field
            header += ' tcool/tff'

        np.savetxt('profiles_{:04d}.txt'.format(int(ds.basename[-4:])), arr,
                   header=header)

# if yt.is_root():
#     last_profile = int(datasets[-1].basename[-4:])
#     with open('last_profile','w') as f:
#         f.write(str(last_profile))
