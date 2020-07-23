import numpy as np
import yt
import scipy.integrate as intg

G = 6.67e-8 * yt.units.cm**3/(yt.units.g * yt.units.s**2)
mh = 1.67e-24
cm_per_kpc = 3.0856e21 * yt.units.cm
gamma = 5/3
mu = 0.59
mu_i = 1.22

Mtot = 2e12 * 1.989e33 # Msun to cgs
Mstar = 2e10 * 1.989e33 * yt.units.g
C = 12
K0 = 2.6 / 6.242e8 # keV cm^2 to erg cm^2
alpha = 0.71

rstar = 3.5 * cm_per_kpc
zstar = 0.325 * cm_per_kpc

# cgs
rho_crit = 1.8788e-29*0.49
r_vir = Rvir = ( 3.0/(4.0*np.pi)*Mtot / (200.*rho_crit) )**(1./3.) * yt.units.cm
Rs = Rvir/C
r_max = Rmax = 2.163*Rs
rho_0 = 200.0*rho_crit * C**3/3.0 / (np.log(1.0+C) - C/(1.0+C)) * yt.units.g/yt.units.cm**3

def M_NFW(r):
    return 4.0*np.pi * rho_0 * Rs**3.0 * (np.log((Rs+r)/Rs) - r/(Rs+r))

def M(r):
    ret = np.where(r <= Rmax, r/Rmax*M_NFW(Rmax), M_NFW(r))
    return ret * yt.units.g

def g(r):
    return G * M(r)/r**2
    
def g_stellar(r, z=0): # assume rcyl = rsph
    return G * Mstar * r * ( r**2 + ( rstar+np.sqrt(z**2+zstar**2) )**2 )**(-3/2)

#def K(r):
#    return K0 * (r/cm_per_kpc)**alpha # erg cm^2
    
#def dPdr(r, P):
#    return -g(r) * 1.2*mu_i*mh*(0.5*P/K(r))**(1/gamma)

# boundary condition on pressure at r_vir
#vcmax2 = G*M(r_max)/r_max
#P_vir = 2*(
#    0.25*mu*mh*vcmax2 / K(r_vir)**(1/gamma)
#     )**( gamma/(gamma-1) ) # dyne cm^-2

#res = intg.solve_ivp(dPdr, (r_vir, 1e-2*cm_per_kpc), (P_vir,)) # ~2-3 times larger than expected

# def point_source_g(r):
#     return G * Mtot / r**2

# def point_source_dPdr(P, r):
#     return -point_source_g(r) * mu_i*mh*(P/K(r))**(1/gamma) 

# # analytic solution to the above
# def point_source_LHS(r):
#     return mu_i*mh*G*Mtot * (cm_per_kpc**alpha/K0)**(1/gamma) / (1+(alpha/gamma)) \
#          * r**( -(1+(alpha/gamma)) )

# def point_source_P(r):
#     const = gamma/(gamma-1) * P_vir**((gamma-1)/gamma) - point_source_LHS(r_vir)
#     print(const)
#     P = point_source_LHS(r) + const
#     P *= (gamma-1)/gamma
#     P = P**( gamma/(gamma-1) )
#     return P




