import numpy as np
import matplotlib.pyplot as plt
import pdb
from pygrackle import \
        chemistry_data, \
        setup_fluid_container
from pygrackle.utilities.physical_constants import \
    mh, \
    gravitational_constant_cgs as G, \
    mass_sun_cgs as Msun, \
    erg_per_keV, \
    kboltz, \
    sec_per_Myr, \
    cm_per_kpc

import h5py
from scipy.interpolate import griddata

# Additional constants
gamma = 5/3
mu = 0.59
mu_i = 1.22

# Parameters
Mtot = 1.1e12 * Msun # cgs
C = 10
K0 = 3.4 * erg_per_keV # keV cm^2 to erg cm^2
alpha = 0.71
n1 = 1.1e-2 # 1/cm^3
n2 = 2.1e-5 # 1/cm^3
zeta1 = 1.2
zeta2 = 2.2
Z = 0.27 # Zsun


# cgs
rho_crit = 1.8788e-29*0.49
r_vir = Rvir = ( 3.0/(4.0*np.pi)*Mtot / (200.*rho_crit) )**(1./3.)
Rs = Rvir/C
r_max = Rmax = 2.163*Rs
rho_0 = 200.0*rho_crit * C**3/3.0 / (np.log(1.0+C) - C/(1.0+C))

def M_NFW(r):
    return 4.0*np.pi * rho_0 * Rs**3.0 * (np.log((Rs+r)/Rs) - r/(Rs+r))

def M(r):
    ret = np.where(r <= Rmax, r/Rmax*M_NFW(Rmax), M_NFW(r))
    return ret

def g(r):
    return G * M(r)/r**2

def K(r_kpc):
    return K0 * (r_kpc)**alpha # erg cm^2

def n_e(r_kpc):
    return 1/np.sqrt( ( n1 * (r_kpc)**-zeta1 )**-2 \
                    + ( n2 * (r_kpc/100)**-zeta2 )**-2  )

def rho(r_kpc):
    return 2*n_e(r_kpc) * mu*mh

def T(r_kpc):
    return K(r_kpc) * n_e(r_kpc)**(2/3) / kboltz # K

# Set solver parameters
chem = chemistry_data()
chem.use_grackle = 1
chem.with_radiative_cooling = 0
chem.primordial_chemistry = 2
chem.metal_cooling = 1
chem.UVbackground = 1
chem.cmb_temperature_floor = 1
chem.grackle_data_file = b"/home/claire/grackle/input/CloudyData_UVB=HM2012.h5"

chem.use_specific_heating_rate = 0
chem.use_volumetric_heating_rate = 0

# Set units
chem.comoving_coordinates = 0 # proper units
chem.a_units = 1.0
chem.a_value = 1.0
# The following are from my Enzo parameter file
chem.density_units = 1.67e-27
chem.length_units = 1638.4 * cm_per_kpc
chem.time_units = sec_per_Myr
chem.velocity_units = chem.a_units \
    * (chem.length_units / chem.a_value) \
    / chem.time_units


# Call convenience function for setting up a fluid container.
# This container holds the solver parameters, units, and fields.
r = np.linspace(30,300,100)

fc = setup_fluid_container(chem,
                           density=rho(r),
                           temperature=T(r),
                           metal_mass_fraction=0.02014*Z,
                           converge=True,
                           tolerance=0.0001)
fc.calculate_temperature()
assert np.allclose(fc['temperature'], T(r))


fc_un = setup_fluid_container(chem,
                              density=rho(r),
                              temperature=T(r),
                              metal_mass_fraction=0.02014*Z,
                              converge=False)

fc2 = setup_fluid_container(chem,
                            density=rho(r),
                            temperature=T(r),
                            metal_mass_fraction=0.02014*Z,
                            converge=False)

#d = h5py.File("equilibrium_table_25_027-Zsun_small2.h5")
d = h5py.File("../equilibrium_table_50_027-Zsun.h5")

dens = np.array(d['indexer']['density'])
temp = np.array(d['indexer']['temperature'])
dd, tt = np.meshgrid(dens, temp)

tab = {}
tab['HI'] = np.array(d['table']['HI'])
tab['HII'] = np.array(d['table']['HII'])
tab['HeI'] = np.array(d['table']['HeI'])
tab['HeII'] = np.array(d['table']['HeII'])
tab['HeIII'] = np.array(d['table']['HeIII'])
tab['H2I'] = np.array(d['table']['H2I'])
tab['H2II'] = np.array(d['table']['H2II'])
tab['HM'] = np.array(d['table']['HM'])
tab['de'] = np.array(d['table']['de'])

d.close()
#import pdb; pdb.set_trace()
label='Table'
for field in fc2.density_fields:
    if field!='density' and field!='metal' and field!='dust':
        fc2[field] = fc2['density'] * griddata((tt.ravel(), dd.ravel()),
                                               tab[field],
                                               (T(r), rho(r)),
                                               method='nearest')

fc.calculate_cooling_time()
fc_un.calculate_cooling_time()
fc2.calculate_cooling_time()

# plt.semilogy(r, rho(r))
# plt.show()
# plt.semilogy(r, T(r))
# plt.show()

tff = np.sqrt( 2 * r*cm_per_kpc / g(r*cm_per_kpc) )
tcool = -fc['cooling_time']*chem.time_units
tcool_un = -fc_un['cooling_time']*chem.time_units
tcool2 = -fc2['cooling_time']*chem.time_units
ratio = tcool / tff
ratio_un = tcool_un / tff
ratio2 = tcool2 / tff

plt.semilogy(r, tcool, label='Equilibrated')
plt.semilogy(r, tcool_un, label='Unequilibrated')
plt.semilogy(r, tcool2, label=label)
plt.scatter(50, 7e16, c='r', label='OoM Estimate')
plt.ylabel(r"$\mathrm{ t_{cool} }$")
plt.xlabel("r [kpc]")
plt.legend()
plt.show()

plt.semilogy(r, ratio, label='Equilibrated')
plt.semilogy(r, ratio_un, label='Unequilibrated')
plt.semilogy(r, ratio2, label=label)
plt.scatter(50, 7e16/7e15, c='r', label='OoM Estimate')
plt.ylabel(r"$\mathrm{ t_{cool}/t_{ff} }$")
plt.xlabel("r [kpc]")
plt.legend()
plt.show()

# I=0
# for field in fc.density_fields :
#     if field != 'density' and field !='dust':     
#         plt.semilogy(r, fc[field], color=f'C{I}', label=field)
#         #plt.semilogy(r, fc_un[field], color=f'C{I}', ls='-', label=field)
#         #plt.semilogy(r, fc2[field], color=f'C{I}', ls=':')
#         #plt.semilogy(r, np.abs((fc2[field]-fc[field])/fc[field]), color=f'C{I}', label=field)
#         I+=1
#         print(field, fc[field][0], fc[field][-1])
# plt.legend(loc='lower right')
# plt.show()

# print(np.array2string(r, separator=','))
# print(np.array2string(fc['HI']/fc['density'], separator=','))
#print(np.array2string(fc['HeII']/fc['density'], separator=','))
# print(np.array2string(rho(r), separator=','))
# print(np.array2string(T(r), separator=','))
