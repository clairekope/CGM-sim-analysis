###################################################
# Compare pressure gradient to NFW profile to check
# hydrostatic equilibrium
###################################################

import yt
import yt.units as u
import numpy as np
import matplotlib.pyplot as plt
from calc_enclosed_mass import *

check_cutout = False


ds = yt.load("DD0000/DD0000")
ds.add_gradient_fields(('gas','pressure')) # creates pressure_gradient_magnitude

dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (27, 'kpc'), (3, 'kpc'))
outskirts = ds.all_data() - dsk
outskirts.set_field_parameter('center', ds.arr([0.5,0.5,0.5],'code_length'))

if check_cutout:
    p = yt.ProjectionPlot(ds, 'x', 'density', width=(200,'kpc'))
    p.set_zlim('density', 5e-5, 2e-2)
    p.save()

    p = yt.ProjectionPlot(ds, 'x', 'density', width=(200,'kpc'),
                          data_source=outskirts)
    p.set_zlim('density', 5e-5, 2e-2)
    p.save('outskirts')

# grad P + rho*g = 0
avg_prof = yt.create_profile(outskirts, 'radius',
                             ['pressure_gradient_magnitude', 'density'],
                             weight_field = 'ones', units={'radius':'kpc'})

r = avg_prof.x
r_bins = avg_prof.x_bins

# functions from save_profiles
Mtot_r = NFW_mass_enclosed(r) +\
         cell_and_particle_mass_enclosed(outskirts, r_bins.v)
g_r = u.G * Mtot_r / np.power(r, 2)
g_r = g_r.to('cm/s**2') + MN_accel(r)

# g_r and grad P are both magnitudes with no sign
HSE_deviation = (avg_prof['pressure_gradient_magnitude'] -\
                 g_r*avg_prof['density']) \
                 / (g_r*avg_prof['density'])

plt.loglog(r, avg_prof['pressure_gradient_magnitude'], label=r'|$\nabla P$|')
plt.loglog(r, g_r*avg_prof['density'], label=r'|$g \rho$|')
plt.axvline(3*206, color='k', ls=':')
plt.legend()
plt.ylabel(r'Quantity (g cm$^{-2}$ s$^{-2}$)')
plt.xlabel('r (kpc)')
plt.savefig('HSE_compare.png')
plt.clf()

plt.semilogx(r, HSE_deviation)
plt.axhline(0, color='k', ls=':')
plt.axvline(3*206, color='k', ls=':')
plt.ylabel(r'(|$\nabla P$| - |$g\rho$|) / |$g\rho$|')
plt.xlabel('r (kpc)')
plt.savefig('HSE_deviation.png')
