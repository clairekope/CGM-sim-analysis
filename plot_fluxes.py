#!/usr/bin/env python
# coding: utf-8

import numpy as np
from astropy.table import Table
from astropy.io.misc.hdf5 import read_table_hdf5
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator


# Fluxes across spherical shells (surfaces, not into volumes)
# Outward is positive.
# 
# Quantities: 
# * thermal_energy
#   - metal thermal_energy
# * thermal energy
# * kinetic energy
#   - radial kinetic energy
#   - tangential kinetic energy
# * potential energy
# * bernoulli energy
# * cooling energy
# * entropy
# 
# Direction:
# * net
# * in
# * out
# 
# Temperatures:
# * [none] is total
# * cold
# * cool
# * warm
# * hot

fid05 = read_table_hdf5("../original_sims/extracted_data/fluxes/fid/fluxes_DD0010_mass_energy.hdf5")
fid15 = read_table_hdf5("../original_sims/extracted_data/fluxes/fid/fluxes_DD0030_mass_energy.hdf5")
fid25 = read_table_hdf5("../original_sims/extracted_data/fluxes/fid/fluxes_DD0050_mass_energy.hdf5")
fid35 = read_table_hdf5("../original_sims/extracted_data/fluxes/fid/fluxes_DD0070_mass_energy.hdf5")

cflow05 = read_table_hdf5("../original_sims/extracted_data/fluxes/cflow/fluxes_DD0010_mass_energy.hdf5")
cflow15 = read_table_hdf5("../original_sims/extracted_data/fluxes/cflow/fluxes_DD0030_mass_energy.hdf5")
cflow25 = read_table_hdf5("../original_sims/extracted_data/fluxes/cflow/fluxes_DD0050_mass_energy.hdf5")
cflow35 = read_table_hdf5("../original_sims/extracted_data/fluxes/cflow/fluxes_DD0070_mass_energy.hdf5")

tctff505 = read_table_hdf5("../original_sims/extracted_data/fluxes/tctff5/fluxes_DD0010_mass_energy.hdf5")
tctff515 = read_table_hdf5("../original_sims/extracted_data/fluxes/tctff5/fluxes_DD0030_mass_energy.hdf5")
tctff525 = read_table_hdf5("../original_sims/extracted_data/fluxes/tctff5/fluxes_DD0050_mass_energy.hdf5")
tctff535 = read_table_hdf5("../original_sims/extracted_data/fluxes/tctff5/fluxes_DD0070_mass_energy.hdf5")

tctff2005 = read_table_hdf5("../original_sims/extracted_data/fluxes/tctff20/fluxes_DD0010_mass_energy.hdf5")
tctff2015 = read_table_hdf5("../original_sims/extracted_data/fluxes/tctff20/fluxes_DD0030_mass_energy.hdf5")
tctff2025 = read_table_hdf5("../original_sims/extracted_data/fluxes/tctff20/fluxes_DD0050_mass_energy.hdf5")
tctff2035 = read_table_hdf5("../original_sims/extracted_data/fluxes/tctff20/fluxes_DD0070_mass_energy.hdf5")

linrot05 = read_table_hdf5("../original_sims/extracted_data/fluxes/linrot/fluxes_DD0010_mass_energy.hdf5")
linrot15 = read_table_hdf5("../original_sims/extracted_data/fluxes/linrot/fluxes_DD0030_mass_energy.hdf5")
linrot25 = read_table_hdf5("../original_sims/extracted_data/fluxes/linrot/fluxes_DD0050_mass_energy.hdf5")
linrot35 = read_table_hdf5("../original_sims/extracted_data/fluxes/linrot/fluxes_DD0070_mass_energy.hdf5")

norot05 = read_table_hdf5("../original_sims/extracted_data/fluxes/norot/fluxes_DD0010_mass_energy.hdf5")
norot15 = read_table_hdf5("../original_sims/extracted_data/fluxes/norot/fluxes_DD0030_mass_energy.hdf5")
norot25 = read_table_hdf5("../original_sims/extracted_data/fluxes/norot/fluxes_DD0050_mass_energy.hdf5")
norot35 = read_table_hdf5("../original_sims/extracted_data/fluxes/norot/fluxes_DD0070_mass_energy.hdf5")


fid_sfh = np.genfromtxt('../original_sims/extracted_data/fid_coarse-bin_sfh.txt')
cflow_sfh = np.genfromtxt('../original_sims/extracted_data/cflow_coarse-bin_sfh.txt')
tctff5_sfh = np.genfromtxt('../original_sims/extracted_data/tctff5_coarse-bin_sfh.txt')
tctff20_sfh = np.genfromtxt('../original_sims/extracted_data/tctff20_coarse-bin_sfh.txt')
linrot_sfh = np.genfromtxt('../original_sims/extracted_data/linrot_coarse-bin_sfh.txt')
norot_sfh = np.genfromtxt('../original_sims/extracted_data/norot_coarse-bin_sfh.txt')


fig, ax = plt.subplots(nrows=5, ncols=4, sharex=True, sharey=True, figsize=(9,10))

ax[0,0].set_title("0.5 Gyr")
ax[0,1].set_title("1.5 Gyr")
ax[0,2].set_title("2.5 Gyr")
ax[0,3].set_title("3.5 Gyr")

ax[0,0].plot(fid05["radius"], fid05["net_cold_mass_flux"]/fid_sfh[10,1], color='C0', label="Fiducial")
ax[0,0].plot(cflow05["radius"], cflow05["net_cold_mass_flux"]/cflow_sfh[10,1], color='C3', ls='-', label="Cooling Flow")
ax[0,0].plot(tctff505["radius"], tctff505["net_cold_mass_flux"]/tctff5_sfh[10,1], color='C2', label=r"$t_{\rm c}/t_{\rm ff}=5$")
ax[0,0].plot(tctff2005["radius"], tctff2005["net_cold_mass_flux"]/tctff20_sfh[10,1], color='C1', ls='-', label=r"$t_{\rm c}/t_{\rm ff}=20$")
ax[0,0].plot(linrot05["radius"], linrot05["net_cold_mass_flux"]/linrot_sfh[10,1], color='C4', ls='-', label="Linear Rotation")
ax[0,0].plot(norot05["radius"], norot05["net_cold_mass_flux"]/norot_sfh[10,1], color='C5', ls='-', label="No Rotation")

ax[0,1].plot(fid15["radius"], fid15["net_cold_mass_flux"]/fid_sfh[30,1], color='C0')
ax[0,1].plot(cflow15["radius"], cflow15["net_cold_mass_flux"]/cflow_sfh[30,1], color='C3', ls='-')
ax[0,1].plot(tctff515["radius"], tctff515["net_cold_mass_flux"]/tctff5_sfh[30,1], color='C2')
ax[0,1].plot(tctff2015["radius"], tctff2015["net_cold_mass_flux"]/tctff20_sfh[30,1], color='C1', ls='-')
ax[0,1].plot(linrot15["radius"], linrot15["net_cold_mass_flux"]/linrot_sfh[30,1], color='C4', ls='-')
ax[0,1].plot(norot15["radius"], norot15["net_cold_mass_flux"]/norot_sfh[30,1], color='C5', ls='-')

ax[0,2].plot(fid25["radius"], fid25["net_cold_mass_flux"]/fid_sfh[50,1], color='C0')
ax[0,2].plot(cflow25["radius"], cflow25["net_cold_mass_flux"]/cflow_sfh[50,1], color='C3', ls='-')
ax[0,2].plot(tctff525["radius"], tctff525["net_cold_mass_flux"]/tctff5_sfh[50,1], color='C2')
ax[0,2].plot(tctff2025["radius"], tctff2025["net_cold_mass_flux"]/tctff20_sfh[50,1], color='C1', ls='-')
ax[0,2].plot(linrot25["radius"], linrot25["net_cold_mass_flux"]/linrot_sfh[50,1], color='C4', ls='-')
ax[0,2].plot(norot25["radius"], norot25["net_cold_mass_flux"]/norot_sfh[50,1], color='C5', ls='-')

ax[0,3].plot(fid35["radius"], fid35["net_cold_mass_flux"]/fid_sfh[70,1], color='C0')
ax[0,3].plot(cflow35["radius"], cflow35["net_cold_mass_flux"]/cflow_sfh[70,1], color='C3', ls='-')
ax[0,3].plot(tctff535["radius"], tctff535["net_cold_mass_flux"]/tctff5_sfh[70,1], color='C2')
ax[0,3].plot(tctff2035["radius"], tctff2035["net_cold_mass_flux"]/tctff20_sfh[70,1], color='C1', ls='-')
ax[0,3].plot(linrot35["radius"], linrot35["net_cold_mass_flux"]/linrot_sfh[70,1], color='C4', ls='-')
ax[0,3].plot(norot35["radius"], norot35["net_cold_mass_flux"]/norot_sfh[70,1], color='C5', ls='-')


ax[1,0].plot(fid05["radius"], fid05["net_cool_mass_flux"]/fid_sfh[10,1], color='C0')
ax[1,0].plot(cflow05["radius"], cflow05["net_cool_mass_flux"]/cflow_sfh[10,1], color='C3', ls='-')
ax[1,0].plot(tctff505["radius"], tctff505["net_cool_mass_flux"]/tctff5_sfh[10,1], color='C2')
ax[1,0].plot(tctff2005["radius"], tctff2005["net_cool_mass_flux"]/tctff20_sfh[10,1], color='C1', ls='-')
ax[1,0].plot(linrot05["radius"], linrot05["net_cool_mass_flux"]/linrot_sfh[10,1], color='C4', ls='-')
ax[1,0].plot(norot05["radius"], norot05["net_cool_mass_flux"]/norot_sfh[10,1], color='C5', ls='-')

ax[1,1].plot(fid15["radius"], fid15["net_cool_mass_flux"]/fid_sfh[30,1], color='C0')
ax[1,1].plot(cflow15["radius"], cflow15["net_cool_mass_flux"]/cflow_sfh[30,1], color='C3', ls='-')
ax[1,1].plot(tctff515["radius"], tctff515["net_cool_mass_flux"]/tctff5_sfh[30,1], color='C2')
ax[1,1].plot(tctff2015["radius"], tctff2015["net_cool_mass_flux"]/tctff20_sfh[30,1], color='C1', ls='-')
ax[1,1].plot(linrot15["radius"], linrot15["net_cool_mass_flux"]/linrot_sfh[30,1], color='C4', ls='-')
ax[1,1].plot(norot15["radius"], norot15["net_cool_mass_flux"]/norot_sfh[30,1], color='C5', ls='-')

ax[1,2].plot(fid25["radius"], fid25["net_cool_mass_flux"]/fid_sfh[50,1], color='C0')
ax[1,2].plot(cflow25["radius"], cflow25["net_cool_mass_flux"]/cflow_sfh[50,1], color='C3', ls='-')
ax[1,2].plot(tctff525["radius"], tctff525["net_cool_mass_flux"]/tctff5_sfh[50,1], color='C2')
ax[1,2].plot(tctff2025["radius"], tctff2025["net_cool_mass_flux"]/tctff20_sfh[50,1], color='C1', ls='-')
ax[1,2].plot(linrot25["radius"], linrot25["net_cool_mass_flux"]/linrot_sfh[50,1], color='C4', ls='-')
ax[1,2].plot(norot25["radius"], norot25["net_cool_mass_flux"]/norot_sfh[50,1], color='C5', ls='-')

ax[1,3].plot(fid35["radius"], fid35["net_cool_mass_flux"]/fid_sfh[70,1], color='C0')
ax[1,3].plot(cflow35["radius"], cflow35["net_cool_mass_flux"]/cflow_sfh[70,1], color='C3', ls='-')
ax[1,3].plot(tctff535["radius"], tctff535["net_cool_mass_flux"]/tctff5_sfh[70,1], color='C2')
ax[1,3].plot(tctff2035["radius"], tctff2035["net_cool_mass_flux"]/tctff20_sfh[70,1], color='C1', ls='-')
ax[1,3].plot(linrot35["radius"], linrot35["net_cool_mass_flux"]/linrot_sfh[70,1], color='C4', ls='-')
ax[1,3].plot(norot35["radius"], norot35["net_cool_mass_flux"]/norot_sfh[70,1], color='C5', ls='-')


ax[2,0].plot(fid05["radius"], fid05["net_warm_mass_flux"]/fid_sfh[10,1], color='C0')
ax[2,0].plot(cflow05["radius"], cflow05["net_warm_mass_flux"]/cflow_sfh[10,1], color='C3', ls='-')
ax[2,0].plot(tctff505["radius"], tctff505["net_warm_mass_flux"]/tctff5_sfh[10,1], color='C2')
ax[2,0].plot(tctff2005["radius"], tctff2005["net_warm_mass_flux"]/tctff20_sfh[10,1], color='C1', ls='-')
ax[2,0].plot(linrot05["radius"], linrot05["net_warm_mass_flux"]/linrot_sfh[10,1], color='C4', ls='-')
ax[2,0].plot(norot05["radius"], norot05["net_warm_mass_flux"]/norot_sfh[10,1], color='C5', ls='-')

ax[2,1].plot(fid15["radius"], fid15["net_warm_mass_flux"]/fid_sfh[30,1], color='C0')
ax[2,1].plot(cflow15["radius"], cflow15["net_warm_mass_flux"]/cflow_sfh[30,1], color='C3', ls='-')
ax[2,1].plot(tctff515["radius"], tctff515["net_warm_mass_flux"]/tctff5_sfh[30,1], color='C2')
ax[2,1].plot(tctff2015["radius"], tctff2015["net_warm_mass_flux"]/tctff20_sfh[30,1], color='C1', ls='-')
ax[2,1].plot(linrot15["radius"], linrot15["net_warm_mass_flux"]/linrot_sfh[30,1], color='C4', ls='-')
ax[2,1].plot(norot15["radius"], norot15["net_warm_mass_flux"]/norot_sfh[30,1], color='C5', ls='-')

ax[2,2].plot(fid25["radius"], fid25["net_warm_mass_flux"]/fid_sfh[50,1], color='C0')
ax[2,2].plot(cflow25["radius"], cflow25["net_warm_mass_flux"]/cflow_sfh[50,1], color='C3', ls='-')
ax[2,2].plot(tctff525["radius"], tctff525["net_warm_mass_flux"]/tctff5_sfh[50,1], color='C2')
ax[2,2].plot(tctff2025["radius"], tctff2025["net_warm_mass_flux"]/tctff20_sfh[50,1], color='C1', ls='-')
ax[2,2].plot(linrot25["radius"], linrot25["net_warm_mass_flux"]/linrot_sfh[50,1], color='C4', ls='-')
ax[2,2].plot(norot25["radius"], norot25["net_warm_mass_flux"]/norot_sfh[50,1], color='C5', ls='-')

ax[2,3].plot(fid35["radius"], fid35["net_warm_mass_flux"]/fid_sfh[70,1], color='C0')
ax[2,3].plot(cflow35["radius"], cflow35["net_warm_mass_flux"]/cflow_sfh[70,1], color='C3', ls='-')
ax[2,3].plot(tctff535["radius"], tctff535["net_warm_mass_flux"]/tctff5_sfh[70,1], color='C2')
ax[2,3].plot(tctff2035["radius"], tctff2035["net_warm_mass_flux"]/tctff20_sfh[70,1], color='C1', ls='-')
ax[2,3].plot(linrot35["radius"], linrot35["net_warm_mass_flux"]/linrot_sfh[70,1], color='C4', ls='-')
ax[2,3].plot(norot35["radius"], norot35["net_warm_mass_flux"]/norot_sfh[70,1], color='C5', ls='-')


ax[3,0].plot(fid05["radius"], fid05["net_hot_mass_flux"]/fid_sfh[10,1], color='C0')
ax[3,0].plot(cflow05["radius"], cflow05["net_hot_mass_flux"]/cflow_sfh[10,1], color='C3', ls='-')
ax[3,0].plot(tctff505["radius"], tctff505["net_hot_mass_flux"]/tctff5_sfh[10,1], color='C2')
ax[3,0].plot(tctff2005["radius"], tctff2005["net_hot_mass_flux"]/tctff20_sfh[10,1], color='C1', ls='-')
ax[3,0].plot(linrot05["radius"], linrot05["net_hot_mass_flux"]/linrot_sfh[10,1], color='C4', ls='-')
ax[3,0].plot(norot05["radius"], norot05["net_hot_mass_flux"]/norot_sfh[10,1], color='C5', ls='-')

ax[3,1].plot(fid15["radius"], fid15["net_hot_mass_flux"]/fid_sfh[30,1], color='C0')
ax[3,1].plot(cflow15["radius"], cflow15["net_hot_mass_flux"]/cflow_sfh[30,1], color='C3', ls='-')
ax[3,1].plot(tctff515["radius"], tctff515["net_hot_mass_flux"]/tctff5_sfh[30,1], color='C2')
ax[3,1].plot(tctff2015["radius"], tctff2015["net_hot_mass_flux"]/tctff20_sfh[30,1], color='C1', ls='-')
ax[3,1].plot(linrot15["radius"], linrot15["net_hot_mass_flux"]/linrot_sfh[30,1], color='C4', ls='-')
ax[3,1].plot(norot15["radius"], norot15["net_hot_mass_flux"]/norot_sfh[30,1], color='C5', ls='-')

ax[3,2].plot(fid25["radius"], fid25["net_hot_mass_flux"]/fid_sfh[50,1], color='C0')
ax[3,2].plot(cflow25["radius"], cflow25["net_hot_mass_flux"]/cflow_sfh[50,1], color='C3', ls='-')
ax[3,2].plot(tctff525["radius"], tctff525["net_hot_mass_flux"]/tctff5_sfh[50,1], color='C2')
ax[3,2].plot(tctff2025["radius"], tctff2025["net_hot_mass_flux"]/tctff20_sfh[50,1], color='C1', ls='-')
ax[3,2].plot(linrot25["radius"], linrot25["net_hot_mass_flux"]/linrot_sfh[50,1], color='C4', ls='-')
ax[3,2].plot(norot25["radius"], norot25["net_hot_mass_flux"]/norot_sfh[50,1], color='C5', ls='-')

ax[3,3].plot(fid35["radius"], fid35["net_hot_mass_flux"]/fid_sfh[70,1], color='C0')
ax[3,3].plot(cflow35["radius"], cflow35["net_hot_mass_flux"]/cflow_sfh[70,1], color='C3', ls='-')
ax[3,3].plot(tctff535["radius"], tctff535["net_hot_mass_flux"]/tctff5_sfh[70,1], color='C2')
ax[3,3].plot(tctff2035["radius"], tctff2035["net_hot_mass_flux"]/tctff20_sfh[70,1], color='C1', ls='-')
ax[3,3].plot(linrot35["radius"], linrot35["net_hot_mass_flux"]/linrot_sfh[70,1], color='C4', ls='-')
ax[3,3].plot(norot35["radius"], norot35["net_hot_mass_flux"]/norot_sfh[70,1], color='C5', ls='-')


ax[4,0].plot(fid05["radius"], fid05["net_mass_flux"]/fid_sfh[10,1], color='C0')
ax[4,0].plot(cflow05["radius"], cflow05["net_mass_flux"]/cflow_sfh[10,1], color='C3', ls='-')
ax[4,0].plot(tctff505["radius"], tctff505["net_mass_flux"]/tctff5_sfh[10,1], color='C2')
ax[4,0].plot(tctff2005["radius"], tctff2005["net_mass_flux"]/tctff20_sfh[10,1], color='C1', ls='-')
ax[4,0].plot(linrot05["radius"], linrot05["net_mass_flux"]/linrot_sfh[10,1], color='C4', ls='-')
ax[4,0].plot(norot05["radius"], norot05["net_mass_flux"]/norot_sfh[10,1], color='C5', ls='-')

ax[4,1].plot(fid15["radius"], fid15["net_mass_flux"]/fid_sfh[30,1], color='C0')
ax[4,1].plot(cflow15["radius"], cflow15["net_mass_flux"]/cflow_sfh[30,1], color='C3', ls='-')
ax[4,1].plot(tctff515["radius"], tctff515["net_mass_flux"]/tctff5_sfh[30,1], color='C2')
ax[4,1].plot(tctff2015["radius"], tctff2015["net_mass_flux"]/tctff20_sfh[30,1], color='C1', ls='-')
ax[4,1].plot(linrot15["radius"], linrot15["net_mass_flux"]/linrot_sfh[30,1], color='C4', ls='-')
ax[4,1].plot(norot15["radius"], norot15["net_mass_flux"]/norot_sfh[30,1], color='C5', ls='-')

ax[4,2].plot(fid25["radius"], fid25["net_mass_flux"]/fid_sfh[50,1], color='C0')
ax[4,2].plot(cflow25["radius"], cflow25["net_mass_flux"]/cflow_sfh[50,1], color='C3', ls='-')
ax[4,2].plot(tctff525["radius"], tctff525["net_mass_flux"]/tctff5_sfh[50,1], color='C2')
ax[4,2].plot(tctff2025["radius"], tctff2025["net_mass_flux"]/tctff20_sfh[50,1], color='C1', ls='-')
ax[4,2].plot(linrot25["radius"], linrot25["net_mass_flux"]/linrot_sfh[50,1], color='C4', ls='-')
ax[4,2].plot(norot25["radius"], norot25["net_mass_flux"]/norot_sfh[50,1], color='C5', ls='-')

ax[4,3].plot(fid35["radius"], fid35["net_mass_flux"]/fid_sfh[70,1], color='C0')
ax[4,3].plot(cflow35["radius"], cflow35["net_mass_flux"]/cflow_sfh[70,1], color='C3', ls='-')
ax[4,3].plot(tctff535["radius"], tctff535["net_mass_flux"]/tctff5_sfh[70,1], color='C2')
ax[4,3].plot(tctff2035["radius"], tctff2035["net_mass_flux"]/tctff20_sfh[70,1], color='C1', ls='-')
ax[4,3].plot(linrot35["radius"], linrot35["net_mass_flux"]/linrot_sfh[70,1], color='C4', ls='-')
ax[4,3].plot(norot35["radius"], norot35["net_mass_flux"]/norot_sfh[70,1], color='C5', ls='-')

ax[0,0].set_xlim(0,200)
ax[0,0].set_ylim(-100,1e4)
ax[0,0].set_yscale('symlog', linthresh=1, linscale=1)

fid = Line2D([0], [0], ls='-', c='C0')
cflow = Line2D([0], [0], ls='-', c='C3')
tctff5 = Line2D([0], [0], ls='-', c='C2')
tctff20 = Line2D([0], [0], ls='-', c='C1')
linrot = Line2D([0], [0], ls='-', c='C4')
norot = Line2D([0], [0], ls='-', c='C5')

ax[0,0].legend([fid, cflow], ["Fiducial","Cooling Flow"], loc="upper right", ncol=1)
ax[1,0].legend([tctff5, tctff20],
               [r"$t_{\rm cool}/t_{\rm ff} = 5$",r"$t_{\rm cool}/t_{\rm ff} =20$"],
               loc="upper right", ncol=1)
ax[2,0].legend([linrot, norot], ["Linear Rotation","No Rotation"], loc="upper right", ncol=1)

for j in range(4):
    ax[4,j].set_xlabel("r [kpc]")
    for i in range(5):
        ax[i,j].axhline(0, c='gray', ls=':')
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(10))
        ax[i,j].fill_between(np.arange(200), -np.ones(200), np.ones(200), color="lightgray")
        ax[i,j].grid(axis='y')
        
ax[0,0].set_ylabel("Net Cold Gas\n"+r"Mass Loading Factor")
ax[1,0].set_ylabel("Net Cool Gas\n"+r"Mass Loading Factor")
ax[2,0].set_ylabel("Net Warm Gas\n"+r"Mass Loading Factor")
ax[3,0].set_ylabel("Net Hot Gas\n"+r"Mass Loading Factor")
ax[4,0].set_ylabel("Net Total Gas\n"+r"Mass Loading Factor", fontweight="bold")

fig.tight_layout()
fig.subplots_adjust(wspace=0.15, hspace=0.15)
fig.savefig("../original_sims/figures/fig_fluxes_all.pdf", dpi=300)




ml20 = np.zeros((5,9))
ml50 = np.zeros((5,9))
ml100 = np.zeros((5,9))

for i in range(9):
    fid = read_table_hdf5(f"../original_sims/extracted_data/fluxes/fid/fluxes_DD00{i}0_mass_energy.hdf5")

    ml20[0,i] = fid["net_cold_mass_flux"][7]/fid_sfh[i*10,1]
    ml20[1,i] = fid["net_cool_mass_flux"][7]/fid_sfh[i*10,1]
    ml20[2,i] = fid["net_warm_mass_flux"][7]/fid_sfh[i*10,1]
    ml20[3,i] = fid["net_hot_mass_flux"][7]/fid_sfh[i*10,1]
    ml20[4,i] = fid["net_mass_flux"][7]/fid_sfh[i*10,1]
    
    ml50[0,i] = fid["net_cold_mass_flux"][22]/fid_sfh[i*10,1]
    ml50[1,i] = fid["net_cool_mass_flux"][22]/fid_sfh[i*10,1]
    ml50[2,i] = fid["net_warm_mass_flux"][22]/fid_sfh[i*10,1]
    ml50[3,i] = fid["net_hot_mass_flux"][22]/fid_sfh[i*10,1]
    ml50[4,i] = fid["net_mass_flux"][22]/fid_sfh[i*10,1]
    
    ml100[0,i] = fid["net_cold_mass_flux"][47]/fid_sfh[i*10,1]
    ml100[1,i] = fid["net_cool_mass_flux"][47]/fid_sfh[i*10,1]
    ml100[2,i] = fid["net_warm_mass_flux"][47]/fid_sfh[i*10,1]
    ml100[3,i] = fid["net_hot_mass_flux"][47]/fid_sfh[i*10,1]
    ml100[4,i] = fid["net_mass_flux"][47]/fid_sfh[i*10,1]

times = np.arange(9)*0.5

fig, ax = plt.subplots(nrows=3, sharex=True, sharey=True, figsize=(4,8))

ax[0].text(0.1, 1e3, "20 kpc", fontsize='large', fontweight='bold',
           verticalalignment='bottom')
ax[0].plot(times, ml20[4], color='k', label="Total")
ax[0].plot(times, ml20[0], ls=":", color='C0', label="Cold")
ax[0].plot(times, ml20[1], ls=(0, (3, 1, 1, 1, 1, 1)), color='C2', label="Cool")
ax[0].plot(times, ml20[2], ls="-.", color='C1', label="Warm")
ax[0].plot(times, ml20[3], ls="--", color='C3', label="Hot")

ax[1].text(0.1, 1e3, "50 kpc", fontsize='large', fontweight='bold',
           verticalalignment='bottom')
ax[1].plot(times, ml50[4], color='k', label="Total")
ax[1].plot(times, ml50[0], ls=":", color='C0', label="Cold")
ax[1].plot(times, ml50[1], ls=(0, (3, 1, 1, 1, 1, 1)), color='C2', label="Cool")
ax[1].plot(times, ml50[2], ls="-.", color='C1', label="Warm")
ax[1].plot(times, ml50[3], ls="--", color='C3', label="Hot")
ax[1].legend(ncol=2)

ax[2].text(0.1, 1e3, "100 kpc", fontsize='large', fontweight='bold',
           verticalalignment='bottom')
ax[2].plot(times, ml100[4], color='k', label="Total")
ax[2].plot(times, ml100[0], ls=":", color='C0', label="Cold")
ax[2].plot(times, ml100[1], ls=(0, (3, 1, 1, 1, 1, 1)), color='C2', label="Cool")
ax[2].plot(times, ml100[2], ls="-.", color='C1', label="Warm")
ax[2].plot(times, ml100[3], ls="--", color='C3', label="Hot")

ax[0].set_yscale('symlog', linthresh=1, linscale=1)

ax[0].set_xlim(0,4)
ax[0].set_ylim(-100,1e4)

ax[2].set_xlabel("t [Gyr]")

for i in range(3):
    ax[i].set_ylabel("Mass Loading Factor")
    ax[i].axhline(0, c='gray', ls=':')
    ax[i].xaxis.set_minor_locator(MultipleLocator(0.5))
    ax[i].fill_between(np.arange(5), -np.ones(5), np.ones(5), color="lightgray")
    ax[i].grid(axis='y')

fig.tight_layout()
fig.savefig("../original_sims/figures/fig_fluxes_fid-ev.pdf", dpi=300)
