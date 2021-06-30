##################################################################
# While save_sfh.py takes the final output and saves the full SFH,
# this script goes through each dataset, calculating the mass of
# stars formed in each star particle (important for feedback)
# and finding the furthest radius at which stars formed in the
# last timestep (to find location of star formation over time)
##################################################################

import yt
import numpy as np
from calc_sfr import calc_sfr
from plot_recent_stars import new_stars

@yt.derived_field(name=('io','particle_mass_formed'), units='g',
                  sampling_type='particle')

def _formed_mass(field, data):
    
    # determine the mass of stars formed between this output and the next.
    x1 = (data.ds.current_time - data[('io','creation_time')])/data[('io','dynamical_time')]
    x2 = (data.ds.current_time + data.ds.quan(50,'Myr') - data[('io','creation_time')])/data[('io','dynamical_time')]
    
    # technically, we stop forming stars more than 12 dyn times from particle creation
    # but that's hard to track with such large timesteps.
    # mass formed, however, should be negligible.
    
    m_form = data[('io','particle_initial_mass')] *\
        ( (1+x1)*np.exp(-x1) - (1+x2)*np.exp(-x2) )
    
    # copying this straight from star_maker2.F
    m_form = np.maximum( np.minimum(m_form, data[('io','particle_mass')]), 0*yt.units.g )
    
    return m_form


datasets = yt.load("../sample_data/fid/DD??6?/DD??6?")
storage = {}
for sto, ds in datasets.piter(dynamic=False, storage=storage):

    stars = ('io', 'particle_type') in ds.field_list
    if stars:
        ds.add_particle_filter('new_stars')
        
    ad = ds.all_data()
    ad.set_field_parameter('normal', ds.arr([0,0,1], 'code_length'))
    
    total_mform = ad.quantities.total_quantity(('io', 'particle_mass_formed')).to('Msun')
    
    farthest_star = ad.quantities.extrema(('new_stars','particle_position_cylindrical_radius')).to('kpc')[1]
    
    sto.result_id = int(ds.basename[-4:])
    sto.result = [total_mform, farthest_star]
    
print(storage)