#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FixedLocator, NullFormatter
import seaborn as sns
import unyt as u
import pickle

# Setup

n_theta = 13
theta = np.round(np.linspace(0, 180, n_theta, endpoint=True))

r_edges = np.linspace(2e-1, 206, 51) * u.kpc
r_centers = r_edges[:-1] + np.diff(r_edges)/2

# Calculate t_ff
# adding in stars as a point source is restricted to r > 30 kpc,
# which is the edge of the stellar distribution at the end of the sim
tff_cut = r_centers > 2*u.kpc

time_myr, fid_mgas, fid_mstar = np.genfromtxt("../extracted_data/fid_masses_over_time.txt", 
                                              usecols=(0,1,2), unpack=True)

Mtot = 1e12 * u.Msun
C = 10

rho_crit = 1.8788e-29*0.49 * u.g/u.cm**3
Rvir = ( 3.0/(4.0*np.pi)*Mtot / (200.*rho_crit) )**(1./3.)
Rs = Rvir/C
rho_0 = 200.0*rho_crit * C**3/3.0 / (np.log(1.0+C) - C/(1.0+C))

M_r = 4.0*np.pi * rho_0 * Rs**3.0 * (np.log((Rs+r_centers)/Rs) - r_centers/(Rs+r_centers))

g_NFW = (u.G*M_r/r_centers**2).to("cm/s**2")

def MN_accel(Mstar_add = None):
    rs = 3.5 * u.kpc
    zs = 0.325 * u.kpc
    MStar = 5.8e10 * u.Msun
    
    if Mstar_add is not None:
        MStar += Mstar_add

    thetacol, radrow = np.meshgrid(theta, r_centers[tff_cut])
    r = radrow*np.sin(thetacol) # cyl radius from sph
    z = radrow*np.cos(thetacol)

    accel_r = u.G*MStar*r/np.power(np.power(r,2) + \
              np.power(rs+np.sqrt(np.power(z,2) + zs**2), 2),3/2)

    accel_z = u.G*MStar*z/np.sqrt(np.power(z,2)+ zs**2) / \
              np.power(np.power(r,2) + np.power(rs+np.sqrt(np.power(z,2)+ zs**2),2),3/2) \
              * (rs+np.sqrt(np.power(z,2) + zs**2))

    g = np.sqrt(np.power(accel_r,2) + np.power(accel_z,2)) # no phi accel
    
    return g.to('cm/s**2')

tff = {}
for t in range(time_myr.size):
    label = f"DD{t:04d}"
    tff[label] = {}
    
    g_MN = MN_accel(fid_mstar[t]*u.Msun)
    
    #g_gas = u.G * fid_mgas[t]*u.Msun / np.power(r_centers[tff_cut], 2)
    
    for i in range(n_theta):
        tff[label][i] = np.zeros(r_centers.size) * np.nan
        g_total = g_NFW + g_MN[:,i]# + g_gas
        tff[label][i][tff_cut] = np.sqrt(2*r_centers[tff_cut] / g_total).to("Gyr")


# Fiducial Late Time (3-4 Gyr) Averages

with open("../extracted_data/fid_hedgehog.pkl","rb") as f:
    fid = pickle.load(f)

ent_meds = np.empty((20, n_theta, 50))
ent_p16s = np.empty_like(ent_meds)
ent_p84s = np.empty_like(ent_meds)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        ent_meds[i,t,:] = fid[f"DD{dd:04}"]["entropy"][t]["med"]
        ent_p16s[i,t,:] = fid[f"DD{dd:04}"]["entropy"][t]["p16"]
        ent_p84s[i,t,:] = fid[f"DD{dd:04}"]["entropy"][t]["p84"]
        
ent_med = np.mean(ent_meds, axis=0)
ent_p16 = np.mean(ent_p16s, axis=0)
ent_p84 = np.mean(ent_p84s, axis=0)


tctff_meds = np.empty((20, n_theta, 50))
tctff_p16s = np.empty_like(tctff_meds)
tctff_p84s = np.empty_like(tctff_meds)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        tctff_meds[i,t,:] = fid[f"DD{dd:04}"]["cooling_time"][t]["med"]/tff[f"DD{dd:04}"][t]
        tctff_p16s[i,t,:] = fid[f"DD{dd:04}"]["cooling_time"][t]["p16"]/tff[f"DD{dd:04}"][t]
        tctff_p84s[i,t,:] = fid[f"DD{dd:04}"]["cooling_time"][t]["p84"]/tff[f"DD{dd:04}"][t]
        
tctff_med = np.mean(tctff_meds, axis=0)
tctff_p16 = np.mean(tctff_p16s, axis=0)
tctff_p84 = np.mean(tctff_p84s, axis=0)


temp_meds = np.empty((20, n_theta, 50))
temp_p16s = np.empty_like(temp_meds)
temp_p84s = np.empty_like(temp_meds)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        temp_meds[i,t,:] = fid[f"DD{dd:04}"]["temperature"][t]["med"]
        temp_p16s[i,t,:] = fid[f"DD{dd:04}"]["temperature"][t]["p16"]
        temp_p84s[i,t,:] = fid[f"DD{dd:04}"]["temperature"][t]["p84"]
        
temp_med = np.mean(temp_meds, axis=0)
temp_p16 = np.mean(temp_p16s, axis=0)
temp_p84 = np.mean(temp_p84s, axis=0)


dens_meds = np.empty((20, n_theta, 50))
dens_p16s = np.empty_like(dens_meds)
dens_p84s = np.empty_like(dens_meds)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        dens_meds[i,t,:] = fid[f"DD{dd:04}"]["density"][t]["med"]
        dens_p16s[i,t,:] = fid[f"DD{dd:04}"]["density"][t]["p16"]
        dens_p84s[i,t,:] = fid[f"DD{dd:04}"]["density"][t]["p84"]
        
dens_med = np.mean(dens_meds, axis=0)
dens_p16 = np.mean(dens_p16s, axis=0)
dens_p84 = np.mean(dens_p84s, axis=0)


mass_meds = np.empty((20, n_theta, 50))
mass_p16s = np.empty_like(mass_meds)
mass_p84s = np.empty_like(mass_meds)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        mass_meds[i,t,:] = fid[f"DD{dd:04}"]["cell_mass"][t]["med"]
        mass_p16s[i,t,:] = fid[f"DD{dd:04}"]["cell_mass"][t]["p16"]
        mass_p84s[i,t,:] = fid[f"DD{dd:04}"]["cell_mass"][t]["p84"]
        
mass_med = np.mean(mass_meds, axis=0)
mass_p16 = np.mean(mass_p16s, axis=0)
mass_p84 = np.mean(mass_p84s, axis=0)


pres_meds = np.empty((20, n_theta, 50))
pres_p16s = np.empty_like(pres_meds)
pres_p84s = np.empty_like(pres_meds)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        pres_meds[i,t,:] = fid[f"DD{dd:04}"]["pressure"][t]["med"]
        pres_p16s[i,t,:] = fid[f"DD{dd:04}"]["pressure"][t]["p16"]
        pres_p84s[i,t,:] = fid[f"DD{dd:04}"]["pressure"][t]["p84"]
        
pres_med = np.mean(pres_meds, axis=0)
pres_p16 = np.mean(pres_p16s, axis=0)
pres_p84 = np.mean(pres_p84s, axis=0)


vel_meds = np.empty((20, n_theta, 50))
vel_p16s = np.empty_like(vel_meds)
vel_p84s = np.empty_like(vel_meds)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        vel_meds[i,t,:] = fid[f"DD{dd:04}"]["radial_velocity"][t]["med"]
        vel_p16s[i,t,:] = fid[f"DD{dd:04}"]["radial_velocity"][t]["p16"]
        vel_p84s[i,t,:] = fid[f"DD{dd:04}"]["radial_velocity"][t]["p84"]
        
vel_med = np.mean(vel_meds, axis=0)
vel_p16 = np.mean(vel_p16s, axis=0)
vel_p84 = np.mean(vel_p84s, axis=0)


# Paper Figure
fig, ax = plt.subplots(nrows=4, sharex=True, sharey=False, figsize=(5,10))

cmap = sns.color_palette("crest_r", as_cmap=True)
colors = [cmap(i) for i in np.abs(theta-90)/90]
zorders = np.concatenate((np.arange(1,n_theta//2+1), np.arange(n_theta//2+1,0,-1)))


for i in range(n_theta):
    ax[0].semilogy(r_centers, ent_med[i], color=colors[i], zorder=zorders[i])
    ax[0].fill_between(r_centers, ent_p16[i], ent_p84[i], alpha=0.2, color=colors[i], zorder=zorders[i])
    ax[2].semilogy(r_centers, tctff_med[i], color=colors[i], zorder=zorders[i])
    ax[2].fill_between(r_centers, tctff_p16[i], tctff_p84[i], alpha=0.2, color=colors[i], zorder=zorders[i])
    ax[3].semilogy(r_centers, mass_med[i], color=colors[i], zorder=zorders[i])
    ax[3].fill_between(r_centers, mass_p16[i], mass_p84[i], alpha=0.2, color=colors[i], zorder=zorders[i])
    ax[1].semilogy(r_centers, pres_med[i], color=colors[i], zorder=zorders[i])
    ax[1].fill_between(r_centers, pres_p16[i], pres_p84[i], alpha=0.2, color=colors[i], zorder=zorders[i])

lines = []
for i in range(n_theta//2+1):
    lines.append(Line2D([],[], color=colors[i]))

for i in range(4):
    ax[i].axvline(20, c='k', ls='--')
    ax[i].xaxis.set_minor_locator(MultipleLocator(5))
    ax[i].tick_params(which='both', top=True, right=True)

ax[2].axhline(5, c='k', ls='--')
ax[2].axhline(20, c='k', ls='--')
    
ax[3].set_xlabel('r  [kpc]', fontsize='large')

ax[0].set_ylabel('K  [keV cm$^2$]', fontsize='large')
ax[2].set_ylabel(r"$t_{\rm c}/t_{\rm ff}$", fontsize='large')
ax[3].set_ylabel(r"$M_{\rm cell}\ \ [\rm M_\odot]$", fontsize='large')
ax[1].set_ylabel(r"P  [erg cm$^{-3}$]", fontsize='large')

ax[0].set_ylim(1e0, 1e6)
ax[2].set_ylim(1e-1, 1e5)
ax[3].set_ylim(1e-2, 1e5)
ax[1].set_ylim(1e-18, 1e-13)
ax[0].set_xlim(0,206)

ax[0].yaxis.set_minor_locator(FixedLocator([1e1, 1e3, 1e5]))
ax[0].yaxis.set_minor_formatter(NullFormatter())
ax[1].yaxis.set_minor_locator(FixedLocator([1e-18, 1e-16, 1e-14]))
ax[1].yaxis.set_minor_formatter(NullFormatter())
ax[2].yaxis.set_minor_locator(FixedLocator([1e0, 1e2, 1e4]))
ax[2].yaxis.set_minor_formatter(NullFormatter())
ax[3].yaxis.set_minor_locator(FixedLocator([1e-1, 1e1, 1e3, 1e5]))
ax[3].yaxis.set_minor_formatter(NullFormatter())

fig.tight_layout()
fig.subplots_adjust(bottom = 0.12)
fig.legend(lines, ['$0^\circ/180^\circ$',
                  '$15^\circ/165^\circ$',
                  '$30^\circ/150^\circ$',
                  '$45^\circ/135^\circ$',
                  '$60^\circ/120^\circ$',
                  '$75^\circ/105^\circ$',
                  '$90^\circ$'],
           loc = 'lower center',
           ncol = 4
         )
fig.savefig("../fig_entropy-pressure-tctff-mass_fid.pdf")


# Variant Late Time Averages

# starting with tctff variants
with open("../extracted_data/tctff5_hedgehog.pkl","rb") as f:
    tctff5 = pickle.load(f)

with open("../extracted_data/tctff20_hedgehog.pkl","rb") as f:
 tctff20 = pickle.load(f)


ent_meds_5 = np.empty((20, n_theta, 50))
ent_p16s_5 = np.empty_like(ent_meds_5)
ent_p84s_5 = np.empty_like(ent_meds_5)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        ent_meds_5[i,t,:] = tctff5[f"DD{dd:04}"]["entropy"][t]["med"]
        ent_p16s_5[i,t,:] = tctff5[f"DD{dd:04}"]["entropy"][t]["p16"]
        ent_p84s_5[i,t,:] = tctff5[f"DD{dd:04}"]["entropy"][t]["p84"]
        
ent_med_5 = np.mean(ent_meds_5, axis=0)
ent_p16_5 = np.mean(ent_p16s_5, axis=0)
ent_p84_5 = np.mean(ent_p84s_5, axis=0)


ent_meds_20 = np.empty((20, n_theta, 50))
ent_p16s_20 = np.empty_like(ent_meds_20)
ent_p84s_20 = np.empty_like(ent_meds_20)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        ent_meds_20[i,t,:] = tctff20[f"DD{dd:04}"]["entropy"][t]["med"]
        ent_p16s_20[i,t,:] = tctff20[f"DD{dd:04}"]["entropy"][t]["p16"]
        ent_p84s_20[i,t,:] = tctff20[f"DD{dd:04}"]["entropy"][t]["p84"]
        
ent_med_20 = np.mean(ent_meds_20, axis=0)
ent_p16_20 = np.mean(ent_p16s_20, axis=0)
ent_p84_20 = np.mean(ent_p84s_20, axis=0)


tctff_meds_5 = np.empty((20, n_theta, 50))
tctff_p16s_5 = np.empty_like(tctff_meds_5)
tctff_p84s_5 = np.empty_like(tctff_meds_5)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        tctff_meds_5[i,t,:] = tctff5[f"DD{dd:04}"]["cooling_time"][t]["med"]/tff[f"DD{dd:04}"][t]
        tctff_p16s_5[i,t,:] = tctff5[f"DD{dd:04}"]["cooling_time"][t]["p16"]/tff[f"DD{dd:04}"][t]
        tctff_p84s_5[i,t,:] = tctff5[f"DD{dd:04}"]["cooling_time"][t]["p84"]/tff[f"DD{dd:04}"][t]
        
tctff_med_5 = np.mean(tctff_meds_5, axis=0)
tctff_p16_5 = np.mean(tctff_p16s_5, axis=0)
tctff_p84_5 = np.mean(tctff_p84s_5, axis=0)


tctff_meds_20 = np.empty((20, n_theta, 50))
tctff_p16s_20 = np.empty_like(tctff_meds_20)
tctff_p84s_20 = np.empty_like(tctff_meds_20)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        tctff_meds_20[i,t,:] = tctff20[f"DD{dd:04}"]["cooling_time"][t]["med"]/tff[f"DD{dd:04}"][t]
        tctff_p16s_20[i,t,:] = tctff20[f"DD{dd:04}"]["cooling_time"][t]["p16"]/tff[f"DD{dd:04}"][t]
        tctff_p84s_20[i,t,:] = tctff20[f"DD{dd:04}"]["cooling_time"][t]["p84"]/tff[f"DD{dd:04}"][t]
        
tctff_med_20 = np.mean(tctff_meds_20, axis=0)
tctff_p16_20 = np.mean(tctff_p16s_20, axis=0)
tctff_p84_20 = np.mean(tctff_p84s_20, axis=0)


mass_meds_5 = np.empty((20, n_theta, 50))
mass_p16s_5 = np.empty_like(mass_meds_5)
mass_p84s_5 = np.empty_like(mass_meds_5)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        mass_meds_5[i,t,:] = tctff5[f"DD{dd:04}"]["cell_mass"][t]["med"]
        mass_p16s_5[i,t,:] = tctff5[f"DD{dd:04}"]["cell_mass"][t]["p16"]
        mass_p84s_5[i,t,:] = tctff5[f"DD{dd:04}"]["cell_mass"][t]["p84"]
        
mass_med_5 = np.mean(mass_meds_5, axis=0)
mass_p16_5 = np.mean(mass_p16s_5, axis=0)
mass_p84_5 = np.mean(mass_p84s_5, axis=0)


mass_meds_20 = np.empty((20, n_theta, 50))
mass_p16s_20 = np.empty_like(mass_meds_20)
mass_p84s_20 = np.empty_like(mass_meds_20)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        mass_meds_20[i,t,:] = tctff20[f"DD{dd:04}"]["cell_mass"][t]["med"]
        mass_p16s_20[i,t,:] = tctff20[f"DD{dd:04}"]["cell_mass"][t]["p16"]
        mass_p84s_20[i,t,:] = tctff20[f"DD{dd:04}"]["cell_mass"][t]["p84"]
        
mass_med_20 = np.mean(mass_meds_20, axis=0)
mass_p16_20 = np.mean(mass_p16s_20, axis=0)
mass_p84_20 = np.mean(mass_p84s_20, axis=0)


vel_meds_5 = np.empty((20, n_theta, 50))
vel_p16s_5 = np.empty_like(vel_meds_5)
vel_p84s_5 = np.empty_like(vel_meds_5)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        vel_meds_5[i,t,:] = tctff5[f"DD{dd:04}"]["radial_velocity"][t]["med"]
        vel_p16s_5[i,t,:] = tctff5[f"DD{dd:04}"]["radial_velocity"][t]["p16"]
        vel_p84s_5[i,t,:] = tctff5[f"DD{dd:04}"]["radial_velocity"][t]["p84"]
        
vel_med_5 = np.mean(vel_meds_5, axis=0)
vel_p16_5 = np.mean(vel_p16s_5, axis=0)
vel_p84_5 = np.mean(vel_p84s_5, axis=0)


vel_meds_20 = np.empty((20, n_theta, 50))
vel_p16s_20 = np.empty_like(vel_meds_20)
vel_p84s_20 = np.empty_like(vel_meds_20)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        vel_meds_20[i,t,:] = tctff20[f"DD{dd:04}"]["radial_velocity"][t]["med"]
        vel_p16s_20[i,t,:] = tctff20[f"DD{dd:04}"]["radial_velocity"][t]["p16"]
        vel_p84s_20[i,t,:] = tctff20[f"DD{dd:04}"]["radial_velocity"][t]["p84"]
        
vel_med_20 = np.mean(vel_meds_20, axis=0)
vel_p16_20 = np.mean(vel_p16s_20, axis=0)
vel_p84_20 = np.mean(vel_p84s_20, axis=0)


# continuing with rotation variants
with open("../extracted_data/linrot_hedgehog.pkl","rb") as f:
    linrot = pickle.load(f)

with open("../extracted_data/norot_hedgehog.pkl","rb") as f:
    norot = pickle.load(f)

ent_meds_lin = np.empty((20, n_theta, 50))
ent_p16s_lin = np.empty_like(ent_meds_lin)
ent_p84s_lin = np.empty_like(ent_meds_lin)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        ent_meds_lin[i,t,:] = linrot[f"DD{dd:04}"]["entropy"][t]["med"]
        ent_p16s_lin[i,t,:] = linrot[f"DD{dd:04}"]["entropy"][t]["p16"]
        ent_p84s_lin[i,t,:] = linrot[f"DD{dd:04}"]["entropy"][t]["p84"]
        
ent_med_lin = np.mean(ent_meds_lin, axis=0)
ent_p16_lin = np.mean(ent_p16s_lin, axis=0)
ent_p84_lin = np.mean(ent_p84s_lin, axis=0)


ent_meds_nor = np.empty((20, n_theta, 50))
ent_p16s_nor = np.empty_like(ent_meds_nor)
ent_p84s_nor = np.empty_like(ent_meds_nor)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        ent_meds_nor[i,t,:] = norot[f"DD{dd:04}"]["entropy"][t]["med"]
        ent_p16s_nor[i,t,:] = norot[f"DD{dd:04}"]["entropy"][t]["p16"]
        ent_p84s_nor[i,t,:] = norot[f"DD{dd:04}"]["entropy"][t]["p84"]
        
ent_med_nor = np.mean(ent_meds_nor, axis=0)
ent_p16_nor = np.mean(ent_p16s_nor, axis=0)
ent_p84_nor = np.mean(ent_p84s_nor, axis=0)


tctff_meds_lin = np.empty((20, n_theta, 50))
tctff_p16s_lin = np.empty_like(tctff_meds_lin)
tctff_p84s_lin = np.empty_like(tctff_meds_lin)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        tctff_meds_lin[i,t,:] = linrot[f"DD{dd:04}"]["cooling_time"][t]["med"]/tff[f"DD{dd:04}"][t]
        tctff_p16s_lin[i,t,:] = linrot[f"DD{dd:04}"]["cooling_time"][t]["p16"]/tff[f"DD{dd:04}"][t]
        tctff_p84s_lin[i,t,:] = linrot[f"DD{dd:04}"]["cooling_time"][t]["p84"]/tff[f"DD{dd:04}"][t]
        
tctff_med_lin = np.mean(tctff_meds_lin, axis=0)
tctff_p16_lin = np.mean(tctff_p16s_lin, axis=0)
tctff_p84_lin = np.mean(tctff_p84s_lin, axis=0)


tctff_meds_nor = np.empty((20, n_theta, 50))
tctff_p16s_nor = np.empty_like(tctff_meds_nor)
tctff_p84s_nor = np.empty_like(tctff_meds_nor)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        tctff_meds_nor[i,t,:] = norot[f"DD{dd:04}"]["cooling_time"][t]["med"]/tff[f"DD{dd:04}"][t]
        tctff_p16s_nor[i,t,:] = norot[f"DD{dd:04}"]["cooling_time"][t]["p16"]/tff[f"DD{dd:04}"][t]
        tctff_p84s_nor[i,t,:] = norot[f"DD{dd:04}"]["cooling_time"][t]["p84"]/tff[f"DD{dd:04}"][t]
        
tctff_med_nor = np.mean(tctff_meds_nor, axis=0)
tctff_p16_nor = np.mean(tctff_p16s_nor, axis=0)
tctff_p84_nor = np.mean(tctff_p84s_nor, axis=0)


vel_meds_lin = np.empty((20, n_theta, 50))
vel_p16s_lin = np.empty_like(vel_meds_lin)
vel_p84s_lin = np.empty_like(vel_meds_lin)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        vel_meds_lin[i,t,:] = linrot[f"DD{dd:04}"]["radial_velocity"][t]["med"]
        vel_p16s_lin[i,t,:] = linrot[f"DD{dd:04}"]["radial_velocity"][t]["p16"]
        vel_p84s_lin[i,t,:] = linrot[f"DD{dd:04}"]["radial_velocity"][t]["p84"]
        
vel_med_lin = np.mean(vel_meds_lin, axis=0)
vel_p16_lin = np.mean(vel_p16s_lin, axis=0)
vel_p84_lin = np.mean(vel_p84s_lin, axis=0)


vel_meds_nor = np.empty((20, n_theta, 50))
vel_p16s_nor = np.empty_like(vel_meds_nor)
vel_p84s_nor = np.empty_like(vel_meds_nor)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        vel_meds_nor[i,t,:] = norot[f"DD{dd:04}"]["radial_velocity"][t]["med"]
        vel_p16s_nor[i,t,:] = norot[f"DD{dd:04}"]["radial_velocity"][t]["p16"]
        vel_p84s_nor[i,t,:] = norot[f"DD{dd:04}"]["radial_velocity"][t]["p84"]
        
vel_med_nor = np.mean(vel_meds_nor, axis=0)
vel_p16_nor = np.mean(vel_p16s_nor, axis=0)
vel_p84_nor = np.mean(vel_p84s_nor, axis=0)

# and now the cooling flow variant
with open("../extracted_data/cflow_hedgehog.pkl","rb") as f:
    cflow = pickle.load(f)

ent_meds_cflow = np.empty((20, n_theta, 50))
ent_p16s_cflow = np.empty_like(ent_meds_cflow)
ent_p84s_cflow = np.empty_like(ent_meds_cflow)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        ent_meds_cflow[i,t,:] = cflow[f"DD{dd:04}"]["entropy"][t]["med"]
        ent_p16s_cflow[i,t,:] = cflow[f"DD{dd:04}"]["entropy"][t]["p16"]
        ent_p84s_cflow[i,t,:] = cflow[f"DD{dd:04}"]["entropy"][t]["p84"]
        
ent_med_cflow = np.mean(ent_meds_cflow, axis=0)
ent_p16_cflow = np.mean(ent_p16s_cflow, axis=0)
ent_p84_cflow = np.mean(ent_p84s_cflow, axis=0)


tctff_meds_cflow = np.empty((20, n_theta, 50))
tctff_p16s_cflow = np.empty_like(tctff_meds_cflow)
tctff_p84s_cflow = np.empty_like(tctff_meds_cflow)

for i in range(20):
    dd = i+60
    for t in range(n_theta):
        tctff_meds_cflow[i,t,:] = cflow[f"DD{dd:04}"]["cooling_time"][t]["med"]/tff[f"DD{dd:04}"][t]
        tctff_p16s_cflow[i,t,:] = cflow[f"DD{dd:04}"]["cooling_time"][t]["p16"]/tff[f"DD{dd:04}"][t]
        tctff_p84s_cflow[i,t,:] = cflow[f"DD{dd:04}"]["cooling_time"][t]["p84"]/tff[f"DD{dd:04}"][t]
        
tctff_med_cflow = np.mean(tctff_meds_cflow, axis=0)
tctff_p16_cflow = np.mean(tctff_p16s_cflow, axis=0)
tctff_p84_cflow = np.mean(tctff_p84s_cflow, axis=0)

# Plot all variants!

fig, ax = plt.subplots(nrows=3, ncols=2, sharex=True, figsize=(6,8))

ind = 6
assert theta[ind] == 90.0

ax[0,0].semilogy(r_centers, fid['DD0000']['entropy'][ind]['med'], color='gray', ls=':')
ax[0,0].semilogy(r_centers, tctff5['DD0000']['entropy'][ind]['med'], color='gray', ls=':')
ax[0,0].semilogy(r_centers, tctff20['DD0000']['entropy'][ind]['med'], color='gray', ls=':')

ax[0,0].semilogy(r_centers, ent_med[ind], color='C0', label='Fiducial')
ax[0,0].fill_between(r_centers, ent_p16[ind], ent_p84[ind], alpha=0.2, color='C0')
ax[0,0].semilogy(r_centers, ent_med_5[ind], color='C2', label=r'$t_{\rm c}/t_{\rm ff} = 5$')
ax[0,0].fill_between(r_centers, ent_p16_5[ind], ent_p84_5[ind], alpha=0.2, color='C2')
ax[0,0].semilogy(r_centers, ent_med_20[ind], color='C1', label=r'$t_{\rm c}/t_{\rm ff} = 20$')
ax[0,0].fill_between(r_centers, ent_p16_20[ind], ent_p84_20[ind], alpha=0.2, color='C1')

# ax[0,0].legend()

ax[0,1].semilogy(r_centers, fid['DD0000']['entropy'][ind]['med'], color='gray', ls=':')

ax[0,1].semilogy(r_centers, ent_med[ind], color='C0', label='Fiducial')
ax[0,1].fill_between(r_centers, ent_p16[ind], ent_p84[ind], alpha=0.2, color='C0')
ax[0,1].semilogy(r_centers, ent_med_lin[ind], color='C4', label='Linear Rotation')
ax[0,1].fill_between(r_centers, ent_p16_lin[ind], ent_p84_lin[ind], alpha=0.2, color='C4')
ax[0,1].semilogy(r_centers, ent_med_nor[ind], color='C5', label='No Rotation')
ax[0,1].fill_between(r_centers, ent_p16_nor[ind], ent_p84_nor[ind], alpha=0.2, color='C5')
#ax[0,1].semilogy(r_centers, ent_p16_nor[ind], color='C5', ls=':')

# ax[0,1].legend()

ax[1,0].semilogy(r_centers, tctff_med[ind], color='C0')
ax[1,0].fill_between(r_centers, tctff_p16[ind], tctff_p84[ind], alpha=0.2, color='C0')
ax[1,0].semilogy(r_centers, tctff_med_5[ind], color='C2')
ax[1,0].fill_between(r_centers, tctff_p16_5[ind], tctff_p84_5[ind], alpha=0.2, color='C2')
ax[1,0].semilogy(r_centers, tctff_med_20[ind], color='C1')
ax[1,0].fill_between(r_centers, tctff_p16_20[ind], tctff_p84_20[ind], alpha=0.2, color='C1')

ax[1,1].semilogy(r_centers, tctff_med[ind], color='C0')
ax[1,1].fill_between(r_centers, tctff_p16[ind], tctff_p84[ind], alpha=0.2, color='C0')
ax[1,1].semilogy(r_centers, tctff_med_lin[ind], color='C4')
ax[1,1].fill_between(r_centers, tctff_p16_lin[ind], tctff_p84_lin[ind], alpha=0.2, color='C4')
ax[1,1].semilogy(r_centers, tctff_med_nor[ind], color='C5')
ax[1,1].fill_between(r_centers, tctff_p16_nor[ind], tctff_p84_nor[ind], alpha=0.2, color='C5')
#ax[1,1].semilogy(r_centers, tctff_p16_nor[ind], color='C5', ls=':')

ax[2,0].plot(r_centers, vel_med[ind], color='C0', label='Fiducial')
ax[2,0].fill_between(r_centers, vel_p16[ind], vel_p84[ind], alpha=0.2, color='C0')
ax[2,0].plot(r_centers, vel_med_5[ind], color='C2', label=r'$t_{\rm c}/t_{\rm ff} = 5$')
ax[2,0].fill_between(r_centers, vel_p16_5[ind], vel_p84_5[ind], alpha=0.2, color='C2')
ax[2,0].plot(r_centers, vel_med_20[ind], color='C1', label=r'$t_{\rm c}/t_{\rm ff} = 20$')
ax[2,0].fill_between(r_centers, vel_p16_20[ind], vel_p84_20[ind], alpha=0.2, color='C1')

ax[2,0].legend(loc="lower right")

ax[2,1].plot(r_centers, vel_med[ind], color='C0', label='Fiducial')
ax[2,1].fill_between(r_centers, vel_p16[ind], vel_p84[ind], alpha=0.2, color='C0')
ax[2,1].plot(r_centers, vel_med_lin[ind], color='C4', label='Linear Rotation')
ax[2,1].fill_between(r_centers, vel_p16_lin[ind], vel_p84_lin[ind], alpha=0.2, color='C4')
ax[2,1].plot(r_centers, vel_med_nor[ind], color='C5', label='No Rotation')
ax[2,1].fill_between(r_centers, vel_p16_nor[ind], vel_p84_nor[ind], alpha=0.2, color='C5')
#ax[2,1].plot(r_centers, vel_p16_nor[ind], color='C5', ls=':')

ax[2,1].legend(loc="lower right")

for i in range(2):
    ax[0,i].set_ylim(1e-2, 1e4)
    ax[1,i].set_ylim(1e-1,1e4)
    ax[2,i].set_ylim(-100,100)
    ax[1,i].axhline(5, color='gray', ls='--')
    ax[1,i].axhline(20, color='gray', ls='--')
    ax[2,i].axhline(0, color='gray', ls='--')
    ax[2,i].set_xlabel("r  [kpc]", fontsize='x-large')
    for j in range(3):
        ax[j,i].axvline(20, color='gray', ls='--')
        ax[j,1].tick_params(labelleft=False)
        
ax[0,0].set_ylabel('K  [keV cm$^2$]', fontsize='x-large')
ax[1,0].set_ylabel(r"$t_{\rm c}/t_{\rm ff}$", fontsize='x-large')
ax[2,0].set_ylabel(r"$v_r$  [km s$^{-1}$]", fontsize='x-large')

ax[0,0].set_xlim(0, 206)
fig.tight_layout()
fig.savefig("../fig_CGM90_entropy-tctff-vel.pdf")
