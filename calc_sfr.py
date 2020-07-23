import yt
import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
yt.enable_parallelism()

#ds = yt.load("DD0601/DD0601")
ds = yt.load("DD0583/DD0583")

#for ds in dataset.piter():
#    if int(ds.basename[-4:]) % 10 == 0 and int(ds.basename[-4]) > 0:

dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (24.5,'kpc'), (2.275, 'kpc'))
#try:
masses = dsk['particle_mass'].in_units('Msun')
formation_time = dsk['creation_time'].in_units('yr')
#except yt.utilities.exceptions.YTFieldNotFound:
#    pass
        
time_range = [0, 2e10] # years
n_bins = 10000
hist, bins = np.histogram(formation_time, bins=n_bins, range=time_range,)
inds = np.digitize(formation_time, bins=bins) # what bin does each time fall in?
time = (bins[:-1] + bins[1:])/2 # center of bin

# Sum masses formed at approx the same time; divide by size of time bins
sfr = np.array([masses[inds == j].sum()/(bins[j+1]-bins[j])
                for j in range(len(time))])
sfr[sfr == 0] = np.nan

plt.plot(time/1e6, sfr)
plt.xlabel('Time  [Myr]')
plt.ylabel('SFR  [M$_\odot$ yr$^{-1}$]')
plt.savefig("SM6_sfr.png")

