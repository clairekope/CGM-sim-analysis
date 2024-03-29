#!/usr/bin/env python
# coding: utf-8

#########################################################################
# Plot multipanel figures comparing sim variants (stored locally) @ 3 Gyr
#########################################################################

import yt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, LinearSegmentedColormap
from matplotlib.ticker import FixedLocator, MultipleLocator, NullFormatter, NullLocator
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

plt.rcParams.update({"font.size":18})

ds_fid60 = yt.load("../sample_data/fid/DD0060/DD0060") # 3 Gyr
ds_cfw60 = yt.load("../sample_data/cflow/DD0060/DD0060")
ds_low60 = yt.load("../sample_data/tctff5/DD0060/DD0060")
ds_hih60 = yt.load("../sample_data/tctff20/DD0060/DD0060")
ds_lin60 = yt.load("../sample_data/linrot/DD0060/DD0060")
ds_nor60 = yt.load("../sample_data/norot/DD0060/DD0060")

# code length is the same in all sims
center = yt.YTQuantity(ds_fid60.quan(0.5,'code_length').to('cm'))
thickness = 2*yt.YTQuantity(3.5,'kpc')
width = yt.YTQuantity(600, 'kpc')
extent = (-width/2, width/2, -width/2, width/2)

fields = ['density','entropy','radial_velocity']

def prep_frb(ds):
    rect = ds.region([center, center, center],
                     [center-thickness/2, center-width/2, center-width/2],
                     [center+thickness/2, center+width/2, center+width/2])

    p = yt.ProjectionPlot(ds, 'x', fields, width=width,
                          data_source=rect, weight_field='ones')
                          
    frb = p.data_source.to_frb(width, 512)

    return frb

frb_cfw60 = prep_frb(ds_cfw60)
frb_low60 = prep_frb(ds_low60)
frb_hih60 = prep_frb(ds_hih60)
frb_lin60 = prep_frb(ds_lin60)
frb_nor60 = prep_frb(ds_nor60)

fig = plt.figure(figsize=(16,9))
grid = ImageGrid(fig, 111, nrows_ncols=(len(fields),5),
               axes_pad=0, label_mode='1', share_all=True,
               cbar_mode='edge', cbar_location='right',
               cbar_pad=0)

grid.axes_llc.tick_params(labelleft=False, labelbottom=False)
for ax in grid:
    ax.tick_params(which='both', axis='both', direction='in')
    ax.xaxis.set_major_locator(FixedLocator([-300,-200,-100,0,100,200,300]))
    ax.xaxis.set_minor_locator(MultipleLocator(20))

d_norm = LogNorm(1e-32, 1e-26)
k_norm = LogNorm(1e0, 1e6)
v_norm = SymLogNorm(1, linscale=0.2, base=10, vmin=-3e3, vmax=3e3)

# with bounds of 1e0 to 1e6, mixed cmap is cleanly divided in half
ent_lo_cmap = plt.get_cmap('crest')(np.linspace(0, 1, 86))
ent_hi_cmap = plt.get_cmap('flare_r')(np.linspace(0, 1, 170))
colors = np.vstack((ent_lo_cmap, ent_hi_cmap))
ent_cmap = LinearSegmentedColormap.from_list("crest_flare", colors)

ax = grid.axes_column[0]
bar = AnchoredSizeBar(ax[2].transData, 100, "100 kpc", 8, 
                       label_top=True, color='white', frameon=False,
                       borderpad=1, size_vertical=5,
                       fontproperties={'size':'x-large',
                                       'weight':'bold'})
ax[2].add_artist(bar)

circle = Circle((0,0), 206, transform=ax[1].transData,
                edgecolor='white', fill=False, ls='--')
ax[1].add_artist(circle)

ax[0].text(0.04, 0.88, "CoolFlow", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_cfw60 = ax[0].imshow(np.array(frb_cfw60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
k_cfw60 = ax[1].imshow(np.array(frb_cfw60['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
v_cfw60 = ax[2].imshow(np.array(frb_cfw60['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)

ax = grid.axes_column[1]
ax[0].text(0.04, 0.88, "LowRatio", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_low60 = ax[0].imshow(np.array(frb_low60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
k_low60 = ax[1].imshow(np.array(frb_low60['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
v_low60 = ax[2].imshow(np.array(frb_low60['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)

ax = grid.axes_column[2]
ax[0].text(0.04, 0.88, "HighRatio", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_hih60 = ax[0].imshow(np.array(frb_hih60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
k_hih60 = ax[1].imshow(np.array(frb_hih60['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
v_hih60 = ax[2].imshow(np.array(frb_hih60['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)
                         
ax = grid.axes_column[3]
ax[0].text(0.04, 0.88, "LinRot", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_lin60 = ax[0].imshow(np.array(frb_lin60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
k_lin60 = ax[1].imshow(np.array(frb_lin60['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
v_lin60 = ax[2].imshow(np.array(frb_lin60['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)

ax = grid.axes_column[4]
ax[0].text(0.04, 0.88, "NoRot", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_nor60 = ax[0].imshow(np.array(frb_nor60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
k_nor60 = ax[1].imshow(np.array(frb_nor60['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
v_nor60 = ax[2].imshow(np.array(frb_nor60['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)


d_cb = fig.colorbar(d_cfw60, cax=grid.cbar_axes[0], extend='both')
k_cb = fig.colorbar(k_cfw60, cax=grid.cbar_axes[1], extend='both')
v_cb = fig.colorbar(v_cfw60, cax=grid.cbar_axes[2], extend='both')

d_cb.set_ticks(FixedLocator([1e-31, 1e-29, 1e-27, 1e-25, 1e-32, 1e-30, 1e-28, 1e-26]))

k_cb.set_ticks(FixedLocator([1e0, 1e2, 1e4, 1e6, 1e-1, 1e1, 1e3, 1e5]))

v_cb.set_ticks(FixedLocator([-1e3,-1e2,-1e1,0,1e1,1e2,1e3]))
v_cb.minorticks_off()

d_cb.set_label(r'Density  [g cm$^{-3}$]')
k_cb.set_label(r'Entropy  [keV cm$^2$]')
v_cb.set_label(r'Radial Velocity  [km/s]')

fig.subplots_adjust(left=0.01, right=0.9, bottom=0.01, top=0.99)
fig.savefig("../fig_edge-comp_poster.pdf", transparent=True)

