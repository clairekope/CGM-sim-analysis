# coding: utf-8
import numpy as np
import yt
G = 6.67e-8 * yt.units.cm**3/(yt.units.g * yt.units.s**2)
mh = 1.67e-24
cm_per_kpc = 3.0856e21
gamma = 5/3
mu = 0.59
mu_i = 1.22

Mtot = 1e12 * 1.989e33 # Msun to cgs
C = 10
K0 = 2.6 / 6.242e8 # keV cm^2 to erg cm^2
alpha = 0.71

# cgs
rho_crit = 1.8788e-29*0.49
r_vir = Rvir = ( 3.0/(4.0*np.pi)*Mtot / (200.*rho_crit) )**(1./3.) * yt.units.cm
Rs = Rvir/C
r_max = Rmax = 2.163*Rs
rho_0 = 200.0*rho_crit * C**3/3.0 / (np.log(1.0+C) - C/(1.0+C)) * yt.units.g/yt.units.cm**3

print(Rvir.to('cm'), np.sqrt(G*Mtot*yt.units.g/Rvir).to('cm/s'))

def M_NFW(r):
    return 4.0*np.pi * rho_0 * Rs**3.0 * (np.log((Rs+r)/Rs) - r/(Rs+r))

def M(r):
    ret = np.where(r <= Rmax, r/Rmax*M_NFW(Rmax), M_NFW(r))
    return ret * yt.units.g

def g(r):
    return G * M(r)/r**2
    
def _tff_mod(field, data):                                                 
    return np.sqrt(2*data['radius']/g(data['radius']))
    
def _tratio(field, data):
    return data['cooling_time']/data['freefall_time']
    
yt.add_field(('gas','freefall_time'), function=_tff_mod, sampling_type='cell', units='s')
yt.add_field(('gas','cooling_freefall_ratio'), function=_tratio, sampling_type='cell', units='1')

#ds = yt.load("DD0000/DD0000")

def g_NFW(r):
    return G * M_NFW(r)/r**2

def g_stellar(r):
    z=0*yt.units.cm
    rs = 3.5*3e21*yt.units.cm
    zs = 0.325*3e21*yt.units.cm
    Mstar = 6.e10*1.989e33*yt.units.g
    return G * Mstar * r / ( r**2 + (rs + np.sqrt(z**2 + zs**2))**2 )**(3/2)
