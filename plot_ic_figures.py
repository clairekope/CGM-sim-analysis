#!/usr/bin/env python
import yt
import unyt as u
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FixedLocator, NullFormatter, MultipleLocator
from calc_enclosed_mass import *


ds = yt.load("../sample_data/fid/DD0000/DD0000")
ds2 = yt.load("../sample_data/fid/DD0060/DD0060")

ad = ds.all_data()
ad.center

# Initial Profiles

profiles = yt.create_profile(ad, "radius", 
                             ["density","temperature","entropy","pressure","cooling_time"],
                             units={"radius":"kpc"}, n_bins=128)

dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (30,'kpc'), (1.3, 'kpc'))
dsk_prof = yt.create_profile(dsk, "radius",
                             ["density","temperature","entropy","pressure","cooling_time"],#,"cell_mass"],
                             units={"radius":"kpc"}, override_bins={"radius":profiles.x})

sph = ds.sphere([0.5,0.5,0.5], (30,'kpc'))
inner_cgm = sph - dsk
inner_cgm.set_field_parameter('center', ds.arr([0.5,0.5,0.5], 'code_length'))
inner_cgm_prof = yt.create_profile(inner_cgm, "radius",
                             ["density","temperature","entropy","pressure","cooling_time"],#,"cell_mass"],
                             units={"radius":"kpc"}, override_bins={"radius":profiles.x})

cgm = ad - dsk
cgm.set_field_parameter('center', ds.arr([0.5,0.5,0.5], 'code_length'))


mass_prof = yt.create_profile(ad, "radius", "cell_mass", units={"radius":"kpc","cell_mass":"Msun"},
                              weight_field=None, accumulation=True,
                              override_bins={"radius":profiles.x})
dsk_mass_prof = yt.create_profile(dsk, "radius", "cell_mass", units={"radius":"kpc","cell_mass":"Msun"},
                                  weight_field=None, accumulation=True,
                                  override_bins={"radius":profiles.x})
cgm_mass_prof = yt.create_profile(cgm, "radius", "cell_mass", units={"radius":"kpc","cell_mass":"Msun"},
                                    weight_field=None, accumulation=True,
                                    override_bins={"radius":profiles.x})

Mtot_r = NFW_mass_enclosed(profiles.x) +          cell_and_particle_mass_enclosed(ad, profiles.x_bins.v)
g_r = u.G * Mtot_r / np.power(profiles.x,2)
g_r = g_r.to('cm/s**2') + MN_accel(profiles.x)
tff = np.sqrt(2*profiles.x/g_r)

fig = plt.figure(figsize=(6,7))
ax1 = plt.subplot2grid(shape=(3,4), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((3,4), (0,2), colspan=2)
ax3 = plt.subplot2grid((3,4), (1,0), colspan=2)
ax4 = plt.subplot2grid((3,4), (1,2), colspan=2)
ax5 = plt.subplot2grid((3,4), (2,0), colspan=2)
ax6 = plt.subplot2grid((3,4), (2,2), colspan=2)

axes = [ax1, ax2, ax3, ax4, ax5, ax6]

ax1.loglog(profiles.x, profiles["density"], lw=2)
ax2.loglog(profiles.x, profiles["temperature"], lw=2)
ax3.loglog(profiles.x, profiles["entropy"], lw=2)
ax4.loglog(profiles.x, profiles["pressure"], lw=2)
ax5.loglog(profiles.x, profiles["cooling_time"].to("Gyr"), lw=2)
ax6.loglog(mass_prof.x, mass_prof["cell_mass"], lw=2)

ax1.loglog(dsk_prof.x[:77], dsk_prof["density"][:77],
           color='C0', ls=':')
ax2.loglog(dsk_prof.x[:77], dsk_prof["temperature"][:77],
           color='C0', ls=':')
ax3.loglog(dsk_prof.x[:77], dsk_prof["entropy"][:77],
           color='C0', ls=':')
ax4.loglog(dsk_prof.x[:77], dsk_prof["pressure"][:77],
           color='C0', ls=':')
ax5.loglog(dsk_prof.x[:77], dsk_prof["cooling_time"].to("Gyr")[:77],
           color='C0', ls=':')

ax1.loglog(inner_cgm_prof.x[37:77], inner_cgm_prof["density"][37:77],
           color='C0', ls='--')
ax2.loglog(inner_cgm_prof.x[37:77], inner_cgm_prof["temperature"][37:77],
           color='C0', ls='--')
ax3.loglog(inner_cgm_prof.x[37:77], inner_cgm_prof["entropy"][37:77],
           color='C0', ls='--')
ax4.loglog(inner_cgm_prof.x[37:77], inner_cgm_prof["pressure"][37:77],
           color='C0', ls='--')
ax5.loglog(inner_cgm_prof.x[37:77], inner_cgm_prof["cooling_time"].to("Gyr")[37:77],
           color='C0', ls='--')
ax6.loglog(cgm_mass_prof.x, cgm_mass_prof["cell_mass"],
           color='C0', ls='--')

ax5.loglog(profiles.x, 10*tff.to('Gyr'), color='gray', ls='-.')

ax1.set_ylim(1e-31, 1e-23)
ax2.set_ylim(5e4, 2e6)
ax3.set_ylim(1e-3, 1e3)
ax4.set_ylim(1e-18, 1e-10)
ax5.set_ylim(1e-7, 1e5)
ax6.set_ylim(1e7,1e11)

ax5.set_xlabel(r"$r\ \ $[kpc]", fontsize='large')
ax6.set_xlabel(r"$r\ \ $[kpc]", fontsize='large')

ax1.set_ylabel(r"$\rho\ \ [\mathrm{g\ cm^{-3}}]$", fontsize='large')
ax2.set_ylabel(r"$T\ \ [\mathrm{K}]$", fontsize='large')
ax3.set_ylabel(r"$K\ \ [\mathrm{keV\ cm^2}]$", fontsize='large')
ax4.set_ylabel(r"$P\ \ [\mathrm{erg\ cm^{-3}}]$", fontsize='large')
ax5.set_ylabel(r"$t_{\rm cool}$  [Gyr]", fontsize='large')
ax6.set_ylabel(r"$M_{\rm enc}(r)\ \ [\mathrm{M_\odot}]$", fontsize='large')

for ax in axes:
    ax.set_xlim(1, 2000)
    ax.axvline(206, c='gray', ls='--', zorder=-2)
    ax.xaxis.set_major_locator(FixedLocator([1e0,1e1,1e2,1e3]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.grid(which='both', alpha=0.4)
    
ax1.yaxis.set_minor_locator(FixedLocator([1e-30, 1e-28, 1e-26, 1e-24]))
ax1.yaxis.set_minor_formatter(NullFormatter())
    
ax3.yaxis.set_minor_locator(FixedLocator([1e-2, 1e0, 1e2]))
ax3.yaxis.set_minor_formatter(NullFormatter())

ax4.yaxis.set_minor_locator(FixedLocator([1e-17, 1e-15, 1e-13, 1e-11]))
ax4.yaxis.set_minor_formatter(NullFormatter())

ax5.yaxis.set_minor_locator(FixedLocator([1e-5, 1e-6, 1e-3, 1e-2, 1e0, 1e1, 1e3, 1e4]))
ax5.yaxis.set_minor_formatter(NullFormatter())

fig.tight_layout()
fig.savefig("../fig_ics.pdf")


# AMR demo

width = ds.quan(600, 'kpc')
width2 = ds.quan(100, 'kpc')
x = np.linspace(-width/2, width/2, 512)
x2 = np.linspace(-width2/2, width2/2, 512)

slc_dx_x2 = yt.SlicePlot(ds2, 'x', 'dy', width=width)
slc_dx_z2 = yt.SlicePlot(ds2, 'z', 'dy', width=width2)

dx_3_x = np.array(slc_dx_x2.data_source.to_frb(width, 512)["dy"].to('pc'))
dx_3_z = np.array(slc_dx_z2.data_source.to_frb(width2, 512)["dy"].to('pc'))

proj_d_x2 = yt.ProjectionPlot(ds2, 'x', 'density', width=width)
proj_d_z2 = yt.ProjectionPlot(ds2, 'z', 'density', width=width2)

d_3_x = np.array(proj_d_x2.data_source.to_frb(width, 512)["density"])
d_3_z = np.array(proj_d_z2.data_source.to_frb(width2, 512)["density"])

cmap = plt.get_cmap('OrRd')
colors = [cmap(i/8+1/16) for i in range(8)]

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(5,7))

# for some reason you can't save levels are a var b/c of float messiness
im = ax[0].imshow(d_3_x, origin="lower", extent=(-width/2,width/2,-width/2,width/2), norm=LogNorm(1e-5,1e-1))
c = ax[0].contour(x, x, dx_3_x, levels=[100, 200, 400, 800, 1600, 3200, 6400, 12800], colors=colors)
# manual placement of contour labels goes in descending order
ax[0].clabel(c, fmt="%0.0f pc", manual=[(-200, -140), # 200, -120
                                        (-180, -110),
                                        (-30, -100),
                                        (0, -50)])

ax[1].imshow(d_3_z, origin="lower", extent=(-width2/2,width2/2,-width2/2,width2/2), norm=LogNorm(1e-5,1e-1))
c = ax[1].contour(x2, x2, dx_3_z, levels=[100, 200, 400, 800, 1600, 3200, 6400, 12800], colors=colors)
c.clabel([200], fmt="%0.0f pc")
ax[1].text(-23,-23,"100 pc", fontsize=10, color=colors[0])

ax[0].set_xlabel('y  [kpc]')
ax[0].set_ylabel('z  [kpc]')

ax[1].set_xlabel('x  [kpc]')
ax[1].set_ylabel('y  [kpc]')

for i in range(2):
    ax[i].tick_params(right=True, top=True)
ax[0].xaxis.set_minor_locator(MultipleLocator(20))
ax[1].xaxis.set_minor_locator(MultipleLocator(5))

fig.subplots_adjust(bottom=0.07, top=0.97, left=0.15, right=0.7, hspace=0.2)
cb_ax = fig.add_axes([0.77, 0.1, 0.05, 0.85])
cb = fig.colorbar(im, cax=cb_ax)
cb.set_label("Surface Density  [g cm$^{-2}$]", fontsize='large')

fig.savefig("../fig_amr_onecol.pdf")

