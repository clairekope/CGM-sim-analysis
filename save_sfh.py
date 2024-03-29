# coding: utf-8
################################################
# Save SFR to file using function in calc_sfr.py
################################################

import glob
import yt
import numpy as np
import calc_sfr

latest_output = sorted(glob.glob("DD????/DD????"))[-1]
ds = yt.load(latest_output)

bins = np.arange(-2.5e7, 4e9+5e7, 5e7) # bins centered on data outputs

if all([field in ds.field_list for field in calc_sfr.required_fields]):

    ds.add_field(name=('nbody','particle_initial_mass'),
                 units='g',
                 sampling_type='particle',
                 function=calc_sfr.field_initial_mass)

time, sfr = calc_sfr.calc_sfr(ds.all_data(),bins)

np.savetxt('sfr.txt',np.column_stack((time/1e6,sfr)),header='Time_Myr SFR_Msun/Yr')
