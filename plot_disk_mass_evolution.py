#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FixedLocator, NullFormatter


# Load disk component masses
fid_mass = np.genfromtxt("../original_sims/extracted_data/fid_masses_over_time.txt", names=True)
cflow_mass = np.genfromtxt("../original_sims/extracted_data/cflow_masses_over_time.txt", names=True)
tctff5_mass = np.genfromtxt("../original_sims/extracted_data/tctff5_masses_over_time.txt", names=True)
tctff20_mass = np.genfromtxt("../original_sims/extracted_data/tctff20_masses_over_time.txt", names=True)
linrot_mass = np.genfromtxt("../original_sims/extracted_data/linrot_masses_over_time.txt", names=True)
norot_mass = np.genfromtxt("../original_sims/extracted_data/norot_masses_over_time.txt", names=True)

# Load star formation history
fid_sfh = np.genfromtxt('../original_sims/extracted_data/fid_coarse-bin_sfh.txt')
cflow_sfh = np.genfromtxt('../original_sims/extracted_data/cflow_coarse-bin_sfh.txt')
tctff5_sfh = np.genfromtxt('../original_sims/extracted_data/tctff5_coarse-bin_sfh.txt')
tctff20_sfh = np.genfromtxt('../original_sims/extracted_data/tctff20_coarse-bin_sfh.txt')
linrot_sfh = np.genfromtxt('../original_sims/extracted_data/linrot_coarse-bin_sfh.txt')
norot_sfh = np.genfromtxt('../original_sims/extracted_data/norot_coarse-bin_sfh.txt')


# Plot disk mass growth and SFH

window = 100 # 0.6 Myr per bin

fig, ax = plt.subplots(nrows=2, ncols=3, sharex=True, figsize=(12,5), gridspec_kw={'height_ratios':[1,0.5]})

for i in range(3):
    ax[0,i].axvline(1, c='gray', ls=':')
    ax[1,i].axvline(1, c='gray', ls=':')

    ax[0,i].axvline(2, c='gray', ls=':')
    ax[1,i].axvline(2, c='gray', ls=':')

    ax[0,i].plot(fid_mass['Time_Myr']/1000, (fid_mass['FormedMass_Msun']+fid_mass['DiskGas_Msun'])/1e9, color='C0', ls='-')
    ax[0,i].plot(fid_mass['Time_Myr']/1000, fid_mass['FormedMass_Msun']/1e9, color='C0', ls='--')
    ax[0,i].plot(fid_mass['Time_Myr']/1000, fid_mass['DiskGas_Msun']/1e9, color='C0', ls=':')
    ax[1,i].plot(fid_sfh[:,0]/1000, fid_sfh[:,1], color='C0', lw=1)

ax[0,0].plot(cflow_mass['Time_Myr']/1000, (cflow_mass['FormedMass_Msun']+cflow_mass['DiskGas_Msun'])/1e9, color='C3', ls='-')
ax[0,0].plot(cflow_mass['Time_Myr']/1000, cflow_mass['FormedMass_Msun']/1e9, color='C3', ls='--')
ax[0,0].plot(cflow_mass['Time_Myr']/1000, cflow_mass['DiskGas_Msun']/1e9, color='C3', ls=':')
ax[1,0].plot(cflow_sfh[:,0]/1000, cflow_sfh[:,1], color='C3', lw=1)

ax[0,1].plot(tctff5_mass['Time_Myr']/1000, (tctff5_mass['FormedMass_Msun']+tctff5_mass['DiskGas_Msun'])/1e9, color='C2', ls='-')
ax[0,1].plot(tctff5_mass['Time_Myr']/1000, tctff5_mass['FormedMass_Msun']/1e9, color='C2', ls='--')
ax[0,1].plot(tctff5_mass['Time_Myr']/1000, tctff5_mass['DiskGas_Msun']/1e9, color='C2', ls=':')
ax[1,1].plot(tctff5_sfh[:,0]/1000, tctff5_sfh[:,1], color='C2', lw=1)

ax[0,1].plot(tctff20_mass['Time_Myr']/1000, (tctff20_mass['FormedMass_Msun']+tctff20_mass['DiskGas_Msun'])/1e9, color='C1', ls='-')
ax[0,1].plot(tctff20_mass['Time_Myr']/1000, tctff20_mass['FormedMass_Msun']/1e9, color='C1', ls='--')
ax[0,1].plot(tctff20_mass['Time_Myr']/1000, tctff20_mass['DiskGas_Msun']/1e9, color='C1', ls=':')
ax[1,1].plot(tctff20_sfh[:,0]/1000, tctff20_sfh[:,1], color='C1', lw=1)

ax[0,2].plot(linrot_mass['Time_Myr']/1000, (linrot_mass['FormedMass_Msun']+linrot_mass['DiskGas_Msun'])/1e9, color='C4', ls='-')
ax[0,2].plot(linrot_mass['Time_Myr']/1000, linrot_mass['FormedMass_Msun']/1e9, color='C4', ls='--')
ax[0,2].plot(linrot_mass['Time_Myr']/1000, linrot_mass['DiskGas_Msun']/1e9, color='C4', ls=':')
ax[1,2].plot(linrot_sfh[:,0]/1000, linrot_sfh[:,1], color='C4', lw=1)

ax[0,2].plot(norot_mass['Time_Myr']/1000, (norot_mass['FormedMass_Msun']+norot_mass['DiskGas_Msun'])/1e9, color='C5', ls='-')
ax[0,2].plot(norot_mass['Time_Myr']/1000, norot_mass['FormedMass_Msun']/1e9, color='C5', ls='--')
ax[0,2].plot(norot_mass['Time_Myr']/1000, norot_mass['DiskGas_Msun']/1e9, color='C5', ls=':')
ax[1,2].plot(norot_sfh[:,0]/1000, norot_sfh[:,1], color='C5', lw=1)

# Set up legend
total = Line2D([0], [0], ls='-', c='dimgray')
stars = Line2D([0], [0], ls='--', c='dimgray')
gas = Line2D([0], [0], ls=':', c='dimgray')

fid = Line2D([0], [0], ls='-', c='C0')
cflow = Line2D([0], [0], ls='-', c='C3')
tctff5 = Line2D([0], [0], ls='-', c='C2')
tctff20 = Line2D([0], [0], ls='-', c='C1')
linrot = Line2D([0], [0], ls='-', c='C4')
norot = Line2D([0], [0], ls='-', c='C5')

ls_fig = ax[0,0].legend([total, stars, gas], ["Disk Mass", "Stellar Mass\nFormed", "Gas Mass"], ncol=1, loc="upper left")
ax[0,0].legend([fid, cflow], ["Fiducial","Cooling Flow"], loc="lower center", ncol=2)
ax[0,1].legend([tctff5, tctff20],
               [r"$t_{\rm cool}/t_{\rm ff} = 5$",r"$t_{\rm cool}/t_{\rm ff} =20$"],
               loc="lower center", ncol=2)
ax[0,2].legend([linrot, norot], ["Linear Rotation","No Rotation"], loc="lower center", ncol=2)
ax[0,0].add_artist(ls_fig)

ax[0,0].set_ylabel(r"$M\ /\ 10^{9}\ {\rm M_\odot}$", fontsize='large')
ax[1,0].set_ylabel(r"SFR$\ \ [{\rm M_\odot/yr}]$", fontsize='large')
for i in range(3):
    ax[1,i].set_xlabel("$t\ \ $[Gyr]", fontsize='large')

ax[0,0].set_xlim(0,4)
for i in range(3):
    ax[0,i].set_ylim(0, 12)
    #ax[1,i].set_ylim(0, 10)
    ax[1,i].set_yscale('log')
    ax[1,i].set_ylim(1e-1,2e1)
    #ax[1,i].yaxis.set_minor_locator(FixedLocator([1e-1, 1e1]))
    #ax[1,i].yaxis.set_minor_formatter(NullFormatter())
    ax[1,i].grid(axis='y', which='major')

for i in range(3):
    for j in range(2):
        ax[j,i].tick_params(right=True, top=True)
for i in range(1,3):
    for j in range(2):
        ax[j,i].tick_params(labelleft=False)

fig.tight_layout()
fig.savefig("../original_sims/figures/fig_mass-ev.pdf")


# Plot disk mass growth


fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, figsize=(12,4), squeeze=False)

for i in range(3):
    ax[0,i].axvline(1, c='gray', ls=':')

    ax[0,i].axvline(2, c='gray', ls=':')

    ax[0,i].plot(fid_mass['Time_Myr']/1000, (fid_mass['FormedMass_Msun']+fid_mass['DiskGas_Msun'])/1e9, color='C0', ls='-')
    ax[0,i].plot(fid_mass['Time_Myr']/1000, fid_mass['FormedMass_Msun']/1e9, color='C0', ls='--')
    ax[0,i].plot(fid_mass['Time_Myr']/1000, fid_mass['DiskGas_Msun']/1e9, color='C0', ls=':')

ax[0,0].plot(cflow_mass['Time_Myr']/1000, (cflow_mass['FormedMass_Msun']+cflow_mass['DiskGas_Msun'])/1e9, color='C3', ls='-')
ax[0,0].plot(cflow_mass['Time_Myr']/1000, cflow_mass['FormedMass_Msun']/1e9, color='C3', ls='--')
ax[0,0].plot(cflow_mass['Time_Myr']/1000, cflow_mass['DiskGas_Msun']/1e9, color='C3', ls=':')

ax[0,1].plot(tctff5_mass['Time_Myr']/1000, (tctff5_mass['FormedMass_Msun']+tctff5_mass['DiskGas_Msun'])/1e9, color='C2', ls='-')
ax[0,1].plot(tctff5_mass['Time_Myr']/1000, tctff5_mass['FormedMass_Msun']/1e9, color='C2', ls='--')
ax[0,1].plot(tctff5_mass['Time_Myr']/1000, tctff5_mass['DiskGas_Msun']/1e9, color='C2', ls=':')

ax[0,1].plot(tctff20_mass['Time_Myr']/1000, (tctff20_mass['FormedMass_Msun']+tctff20_mass['DiskGas_Msun'])/1e9, color='C1', ls='-')
ax[0,1].plot(tctff20_mass['Time_Myr']/1000, tctff20_mass['FormedMass_Msun']/1e9, color='C1', ls='--')
ax[0,1].plot(tctff20_mass['Time_Myr']/1000, tctff20_mass['DiskGas_Msun']/1e9, color='C1', ls=':')

ax[0,2].plot(linrot_mass['Time_Myr']/1000, (linrot_mass['FormedMass_Msun']+linrot_mass['DiskGas_Msun'])/1e9, color='C4', ls='-')
ax[0,2].plot(linrot_mass['Time_Myr']/1000, linrot_mass['FormedMass_Msun']/1e9, color='C4', ls='--')
ax[0,2].plot(linrot_mass['Time_Myr']/1000, linrot_mass['DiskGas_Msun']/1e9, color='C4', ls=':')

ax[0,2].plot(norot_mass['Time_Myr']/1000, (norot_mass['FormedMass_Msun']+norot_mass['DiskGas_Msun'])/1e9, color='C5', ls='-')
ax[0,2].plot(norot_mass['Time_Myr']/1000, norot_mass['FormedMass_Msun']/1e9, color='C5', ls='--')
ax[0,2].plot(norot_mass['Time_Myr']/1000, norot_mass['DiskGas_Msun']/1e9, color='C5', ls=':')

# Set up legend
total = Line2D([0], [0], ls='-', c='dimgray')
stars = Line2D([0], [0], ls='--', c='dimgray')
gas = Line2D([0], [0], ls=':', c='dimgray')

fid = Line2D([0], [0], ls='-', c='C0')
cflow = Line2D([0], [0], ls='-', c='C3')
tctff5 = Line2D([0], [0], ls='-', c='C2')
tctff20 = Line2D([0], [0], ls='-', c='C1')
linrot = Line2D([0], [0], ls='-', c='C4')
norot = Line2D([0], [0], ls='-', c='C5')

ls_fig = ax[0,0].legend([total, stars, gas], ["Disk Mass", "Stellar Mass\nFormed", "Gas Mass"], ncol=1, loc="upper left")
ax[0,0].legend([fid, cflow], ["Fiducial","Cooling Flow"], loc="lower center", ncol=2)
ax[0,1].legend([tctff5, tctff20],
               [r"$t_{\rm cool}/t_{\rm ff} = 5$",r"$t_{\rm cool}/t_{\rm ff} =20$"],
               loc="lower center", ncol=2)
ax[0,2].legend([linrot, norot], ["Linear Rotation","No Rotation"], loc="lower center", ncol=2)
ax[0,0].add_artist(ls_fig)

ax[0,0].set_ylabel(r"$M\ /\ 10^{9}\ {\rm M_\odot}$", fontsize='large')
for i in range(3):
    ax[0,i].set_xlabel("$t\ \ $[Gyr]", fontsize='large')

ax[0,0].set_xlim(0,4)
for i in range(3):
    ax[0,i].set_ylim(0, 12)
    #ax[1,i].set_ylim(0, 10)
    #ax[1,i].yaxis.set_minor_locator(FixedLocator([1e-1, 1e1]))
    #ax[1,i].yaxis.set_minor_formatter(NullFormatter())

for i in range(3):
        ax[0,i].tick_params(right=True, top=True)
for i in range(1,3):
        ax[0,i].tick_params(labelleft=False)

fig.tight_layout()
fig.savefig("../original_sims/figures/fig_mass-ev_noSFR.png")


# Late-Time Change in disk mass

ref = 40
end = 81
fid_ref = fid_mass[ref]
cflow_ref = cflow_mass[ref]
tctff5_ref = tctff5_mass[ref]
tctff20_ref = tctff20_mass[ref]
linrot_ref = linrot_mass[ref]
norot_ref = norot_mass[ref]

fid_late_time = fid_mass["Time_Myr"][ref:end]/1000
fid_late_gas = fid_mass["DiskGas_Msun"][ref:end]-fid_ref["DiskGas_Msun"]
fid_late_form = fid_mass["FormedMass_Msun"][ref:end]-fid_ref["FormedMass_Msun"]
fid_late_sfr = fid_sfh[(fid_sfh[:,0]/1000 >= 2) & (fid_sfh[:,0]/1000 <= 4)]

cflow_late_time = cflow_mass["Time_Myr"][ref:end]/1000
cflow_late_gas = cflow_mass["DiskGas_Msun"][ref:end]-cflow_ref["DiskGas_Msun"]
cflow_late_form = cflow_mass["FormedMass_Msun"][ref:end]-cflow_ref["FormedMass_Msun"]
cflow_late_sfr = cflow_sfh[(cflow_sfh[:,0]/1000 >= 2) & (cflow_sfh[:,0]/1000 <= 4)]

tctff5_late_time = tctff5_mass["Time_Myr"][ref:end]/1000
tctff5_late_gas = tctff5_mass["DiskGas_Msun"][ref:end]-tctff5_ref["DiskGas_Msun"]
tctff5_late_form = tctff5_mass["FormedMass_Msun"][ref:end]-tctff5_ref["FormedMass_Msun"]
tctff5_late_sfr = tctff5_sfh[(tctff5_sfh[:,0]/1000 >= 2) & (tctff5_sfh[:,0]/1000 <= 4)]

tctff20_late_time = tctff20_mass["Time_Myr"][ref:end]/1000
tctff20_late_gas = tctff20_mass["DiskGas_Msun"][ref:end]-tctff20_ref["DiskGas_Msun"]
tctff20_late_form = tctff20_mass["FormedMass_Msun"][ref:end]-tctff20_ref["FormedMass_Msun"]
tctff20_late_sfr = tctff20_sfh[(tctff20_sfh[:,0]/1000 >= 2) & (tctff20_sfh[:,0]/1000 <= 4)]

linrot_late_time = linrot_mass["Time_Myr"][ref:end]/1000
linrot_late_gas = linrot_mass["DiskGas_Msun"][ref:end]-linrot_ref["DiskGas_Msun"]
linrot_late_form = linrot_mass["FormedMass_Msun"][ref:end]-linrot_ref["FormedMass_Msun"]
linrot_late_sfr = linrot_sfh[(linrot_sfh[:,0]/1000 >= 2) & (linrot_sfh[:,0]/1000 <= 4)]

norot_late_time = norot_mass["Time_Myr"][ref:end]/1000
norot_late_gas = norot_mass["DiskGas_Msun"][ref:end]-norot_ref["DiskGas_Msun"]
norot_late_form = norot_mass["FormedMass_Msun"][ref:end]-norot_ref["FormedMass_Msun"]
norot_late_sfr = norot_sfh[(norot_sfh[:,0]/1000 >= 2) & (norot_sfh[:,0]/1000 <= 4)]


# Plot late-time growth

fig, ax = plt.subplots(nrows=2, sharex=True, figsize=(5,6))

ax[0].plot(fid_late_time,
         fid_late_form/1e8,
         c="C0", label="Fiducial")
ax[0].plot(cflow_late_time,
         cflow_late_form/1e8,
         c='C3', label="Cooling Flow")
ax[0].plot(tctff5_late_time,
         tctff5_late_form/1e8,
         c="C2", label=r"$t_{\rm cool}/t_{\rm ff} = 5$")
ax[0].plot(tctff20_late_time,
         tctff20_late_form/1e8,
         c="C1")
ax[0].plot(linrot_late_time,
         linrot_late_form/1e8,
         c="C4")
ax[0].plot(norot_late_time,
         norot_late_form/1e8,
         c="C5")

ax[1].plot(fid_late_time,
         fid_late_gas/1e8,
         c="C0")
ax[1].plot(cflow_late_time,
         cflow_late_gas/1e8,
         c='C3')
ax[1].plot(tctff5_late_time,
         tctff5_late_gas/1e8,
         c="C2")
ax[1].plot(tctff20_late_time,
         tctff20_late_gas/1e8,
         c="C1", label=r"$t_{\rm cool}/t_{\rm ff} = 20$")
ax[1].plot(linrot_late_time,
         linrot_late_gas/1e8,
         c="C4", label="Linear Rotation")
ax[1].plot(norot_late_time,
         norot_late_gas/1e8,
         c="C5", label="No Rotation")

ax[1].axhline(0, c='gray', ls='--')

ax[0].set_xlim(2, 4)

ax[0].set_ylim(0, 18)
ax[1].set_ylim(-8.5, 6)

ax[1].set_xlabel("$t\ \ $[Gyr]", fontsize='large')

ax[0].set_ylabel(r"Stellar Mass ($\rm 10^8\ M_\odot$)", fontsize='large')
ax[1].set_ylabel(r"Gas Mass ($\rm 10^8\ M_\odot$)", fontsize='large')

ax[0].legend(loc='upper left')
ax[1].legend(loc='lower left')

for i in range(2):
    ax[i].tick_params(which='both', right=True, top=True)

ax[0].yaxis.set_major_locator(MultipleLocator(2))

fig.tight_layout()
fig.savefig("../original_sims/figures/fig_mass-ev_late.pdf")


# Calculate late-time rates

fid_form_rate = np.gradient(fid_late_form, 5e7)
fid_gas_rate = np.gradient(fid_late_gas, 5e7) + fid_form_rate

cflow_form_rate = np.gradient(cflow_late_form, 5e7)
cflow_gas_rate = np.gradient(cflow_late_gas, 5e7) + cflow_form_rate

tctff5_form_rate = np.gradient(tctff5_late_form, 5e7)
tctff5_gas_rate = np.gradient(tctff5_late_gas, 5e7) + tctff5_form_rate

tctff20_form_rate = np.gradient(tctff20_late_form, 5e7)
tctff20_gas_rate = np.gradient(tctff20_late_gas, 5e7) + tctff20_form_rate

linrot_form_rate = np.gradient(linrot_late_form, 5e7)
linrot_gas_rate = np.gradient(linrot_late_gas, 5e7) + linrot_form_rate

norot_form_rate = np.gradient(norot_late_form, 5e7)
norot_gas_rate = np.gradient(norot_late_gas, 5e7) + norot_form_rate


# Plot rates

fig, ax = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(12,5))

ax[0,0].plot(fid_late_sfr[:,0]/1000, fid_late_sfr[:,1], color='dimgray')
ax[0,0].plot(fid_late_time, fid_gas_rate, lw=3)
ax[0,0].hlines(fid_gas_rate[ 0:20].mean(), 2, 3, color='C0', lw=1)
ax[0,0].hlines(fid_gas_rate[20:41].mean(), 3, 4, color='C0', lw=1)
ax[0,0].axhline(0, color='k', ls='--')
ax[0,0].set_ylabel(r"$\Delta M$  [M$_\odot$ yr$^{-1}$]", labelpad=1.2, fontsize='large')
ax[0,0].set_title("Fiducial ", loc='right', y=0.85)

ax[1,0].plot(cflow_late_sfr[:,0]/1000, cflow_late_sfr[:,1], color='dimgray')
ax[1,0].plot(cflow_late_time, cflow_gas_rate, lw=3, c='C3')
ax[1,0].hlines(cflow_gas_rate[ 0:20].mean(), 2, 3, color='C3', lw=1)
ax[1,0].hlines(cflow_gas_rate[20:41].mean(), 3, 4, color='C3', lw=1)
ax[1,0].axhline(0, color='k', ls='--')
ax[1,0].set_ylabel(r"$\Delta M$  [M$_\odot$ yr$^{-1}$]", labelpad=1.2, fontsize='large')
ax[1,0].set_title("Cooling Flow ", loc='right', y=0.85)

ax[0,1].plot(tctff5_late_sfr[:,0]/1000, tctff5_late_sfr[:,1], color='dimgray')
ax[0,1].plot(tctff5_late_time, tctff5_gas_rate, lw=3, c='C2')
ax[0,1].hlines(tctff5_gas_rate[ 0:20].mean(), 2, 3, color='C2', lw=1)
ax[0,1].hlines(tctff5_gas_rate[20:41].mean(), 3, 4, color='C2', lw=1)
ax[0,1].axhline(0, color='k', ls='--')
ax[0,1].set_title(r"$t_{\rm c}/t_{\rm ff} = 5\ $ ", loc='right', y=0.85)

ax[1,1].plot(tctff20_late_sfr[:,0]/1000, tctff20_late_sfr[:,1], color='dimgray')
ax[1,1].plot(tctff20_late_time, tctff20_gas_rate, lw=3, c='C1')
ax[1,1].hlines(tctff20_gas_rate[ 0:20].mean(), 2, 3, color='C1', lw=1)
ax[1,1].hlines(tctff20_gas_rate[20:41].mean(), 3, 4, color='C1', lw=1)
ax[1,1].axhline(0, color='k', ls='--')
ax[1,1].set_title(r"$t_{\rm c}/t_{\rm ff} = 20\ $ ", loc='right', y=0.85)

ax[0,2].plot(linrot_late_sfr[:,0]/1000, linrot_late_sfr[:,1], color='dimgray')
ax[0,2].plot(linrot_late_time, linrot_gas_rate, lw=3, c='C4')
ax[0,2].hlines(linrot_gas_rate[ 0:20].mean(), 2, 3, color='C4', lw=1)
ax[0,2].hlines(linrot_gas_rate[20:41].mean(), 3, 4, color='C4', lw=1)
ax[0,2].axhline(0, color='k', ls='--')
ax[0,2].set_title("Linear Rotation ", loc='right', y=0.85)

ax[1,2].plot(norot_late_sfr[:,0]/1000, norot_late_sfr[:,1], color='dimgray')
ax[1,2].plot(norot_late_time, norot_gas_rate, lw=3, c='C5')
ax[1,2].hlines(norot_gas_rate[ 0:20].mean(), 2, 3, color='C5', lw=1)
ax[1,2].hlines(norot_gas_rate[20:41].mean(), 3, 4, color='C5', lw=1)
ax[1,2].axhline(0, color='k', ls='--')
ax[1,2].set_title("No Rotation ", loc='right', y=0.85)

ax[0,0].set_xlim(2, 4)
ax[0,0].set_ylim(-1.4,2.4)

SFR = Line2D([0], [0], ls='-', c='dimgray')
mdisk = Line2D([0], [0], ls='-', lw=3, c='C0')

ax[0,0].legend([SFR,mdisk], ['SFR', r'$M_{\rm disk}$'],
               ncol=2, loc="lower center")

for i in range(3):
    ax[1,i].set_xlabel('$t$  [Gyr]', fontsize='large')

fig.tight_layout()
fig.savefig("../original_sims/figures/fig_net-gas-mass.pdf", dpi=300)
