#####################################################
# Calculate the star formation rate in provided bins
# When run standalone, will save as png (but not txt)
#####################################################

import yt
from yt.fields.api import ValidateDataField
import numpy as np

def field_initial_mass(field, data):
    # ripped from star_maker2.F, assuming StarMassEjectionFraction = 0.25
    # Also assuming all particles are stars!
    
    time_frac = (data.ds.current_time - data[('nbody','creation_time')]) \
                / data[('nbody','dynamical_time')]

    return data[('nbody','particle_mass')] \
           / (1 - 0.25*(1 - (1+time_frac)*np.exp(-time_frac)))

def field_formed_mass(field, data):
    # After back-calculating the initial particle mass, we can find the
    # mass of stars that have formed as used in the FB algorithm

    time_frac = (data.ds.current_time - data[('nbody','creation_time')]) \
                  / data[('nbody','dynamical_time')]
    return data[('nbody','particle_initial_mass')] * time_frac

def calc_sfr(obj, year_bins):

    masses = obj[('nbody','particle_initial_mass')].in_units('Msun')
    formation_time = obj[('nbody','creation_time')].in_units('yr')
        
    inds = np.digitize(formation_time, bins=year_bins) # what bin does each time fall in?
    time = (year_bins[:-1] + year_bins[1:])/2 # center of bin
    
    # Sum masses formed at approx the same time; divide by size of time bins
    sfr = np.array([masses[inds == j].sum()/(year_bins[j+1]-year_bins[j])
                    for j in range(len(time))])
    sfr[sfr == 0] = np.nan

    return time, sfr


yt.add_field(
    ('nbody','particle_initial_mass'),
    function = field_initial_mass,
    units='g',
    sampling_type='particle',
    validators=[ValidateDataField(('nbody','creation_time')),
                ValidateDataField(("nbody","dynamical_time")),
                ValidateDataField(("nbody","particle_mass"))]
)

yt.add_field(
    ('nbody','particle_formed_mass'),
    function = field_formed_mass,
    units='g',
    sampling_type='particle',
    validators=[ValidateDataField(('nbody','creation_time')),
                ValidateDataField(("nbody","dynamical_time")),
                ValidateDataField(("nbody","particle_mass"))]
)


if __name__=="__main__":

    import glob
    import os
    import matplotlib; matplotlib.use('agg')
    import matplotlib.pyplot as plt
    
    latest_output = sorted(glob.glob("DD????/DD????"))[-1]

    ds = yt.load(latest_output)

    dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (24.5,'kpc'), (2.275, 'kpc'))

    bins = ds.arr(np.linspace(0, 2e10, 10001), "yr")
    time, sfr = calc_sfr(dsk, bins) 

    cwd = os.getcwd()

    plt.plot(time.to("Myr"), sfr)
    plt.ylim(0, 10)
    plt.xlabel('Time  [Myr]')
    plt.ylabel('SFR  [M$_\odot$ yr$^{-1}$]')
    plt.title(os.path.basename(cwd))
    plt.savefig("sfr.png")

