# coding: utf-8
################################################
# Save SFR to file using function in calc_sfr.py
################################################

import glob
import yt
import numpy as np
from calc_sfr import calc_sfr

latest_output = sorted(glob.glob("DD????/DD????"))[-1]
ds = yt.load(latest_output)

bins = np.linspace(0, 6e9, 10000)
time, sfr = calc_sfr(ds.all_data(),bins)

np.savetxt('sfr.txt',np.column_stack((time/1e6,sfr)),header='Time_Myr SFR_Msun/Yr')
