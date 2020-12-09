import yt
import glob
import os
import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
yt.enable_parallelism()

def _initial_mass(field, data):
    # ripped from star_maker2.F, assuming StarMassEjectionFraction = 0.25
    time_frac = (ds.current_time - data[('io','creation_time')]) \
                / data[('io','dynamical_time')]
    return data[('io','particle_mass')] \
        / (1 - 0.25*(1 - (1+time_frac)*np.exp(-time_frac)))

yt.add_field('particle_initial_mass', function=_initial_mass, units='g')

def calc_sfr(obj, year_bounds, n_bins):

    masses = obj['particle_initial_mass'].in_units('Msun')
    formation_time = obj['creation_time'].in_units('yr')
        
    hist, bins = np.histogram(formation_time, bins=n_bins, range=year_bounds)
    inds = np.digitize(formation_time, bins=bins) # what bin does each time fall in?
    time = (bins[:-1] + bins[1:])/2 # center of bin
    
    # Sum masses formed at approx the same time; divide by size of time bins
    sfr = np.array([masses[inds == j].sum()/(bins[j+1]-bins[j])
                    for j in range(len(time))])
    sfr[sfr == 0] = np.nan

    return time, sfr

if __name__=="__main__":

    latest_output = sorted(glob.glob("DD????/DD????"))[-1]

    ds = yt.load(latest_output)

    dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (24.5,'kpc'), (2.275, 'kpc'))

    time, sfr = calc_sfr(dsk, (0, 2e10), 10000) 

    cwd = os.getcwd()

    plt.plot(time/1e6, sfr)
    plt.ylim(0, 10)
    plt.xlabel('Time  [Myr]')
    plt.ylabel('SFR  [M$_\odot$ yr$^{-1}$]')
    plt.title(os.path.basename(cwd))
    plt.savefig("sfr.png")

