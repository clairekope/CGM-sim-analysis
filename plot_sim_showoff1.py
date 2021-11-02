#!/usr/bin/env python
# coding: utf-8

###############################################################
# Plot multipanel figures showing evolution of the Fiducial sim
###############################################################

import yt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, LinearSegmentedColormap
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

ds_fid20 = yt.load("DD0020/DD0020") # 1 Gyr
ds_fid40 = yt.load("DD0040/DD0040") # 2 Gyr
ds_fid60 = yt.load("DD0060/DD0060") # 3 Gyr
ds_fid80 = yt.load("DD0080/DD0080") # 4 Gyr

# code length is the same in all sims
center = yt.YTQuantity(ds_fid20.quan(0.5,'code_length').to('cm'))
thickness = 2*yt.YTQuantity(3.5,'kpc')
width = yt.YTQuantity(600, 'kpc')
extent = (-width/2, width/2, -width/2, width/2)

fields = ['density','temperature','entropy','pressure','radial_velocity']

rect_fid20 = ds_fid20.region([center, center, center],
                             [center-thickness/2, center-width/2, center-width/2],
                             [center+thickness/2, center+width/2, center+width/2])

rect_fid40 = ds_fid40.region([center, center, center],
                             [center-thickness/2, center-width/2, center-width/2],
                             [center+thickness/2, center+width/2, center+width/2])

rect_fid60 = ds_fid60.region([center, center, center],
                             [center-thickness/2, center-width/2, center-width/2],
                             [center+thickness/2, center+width/2, center+width/2])

rect_fid80 = ds_fid80.region([center, center, center],
                             [center-thickness/2, center-width/2, center-width/2],
                             [center+thickness/2, center+width/2, center+width/2])


p_fid20 = yt.ProjectionPlot(ds_fid20, 'x', fields, width=width,
                            data_source=rect_fid20, weight_field='ones')
frb_fid20 = p_fid20.data_source.to_frb(width, 512)

p_fid40 = yt.ProjectionPlot(ds_fid40, 'x', fields, width=width,
                            data_source=rect_fid40, weight_field='ones')
frb_fid40 = p_fid40.data_source.to_frb(width, 512)

p_fid60 = yt.ProjectionPlot(ds_fid60, 'x', fields, width=width,
                            data_source=rect_fid60, weight_field='ones')
frb_fid60 = p_fid60.data_source.to_frb(width, 512)

p_fid80 = yt.ProjectionPlot(ds_fid80, 'x', fields, width=width,
                            data_source=rect_fid80, weight_field='ones')
frb_fid80 = p_fid80.data_source.to_frb(width, 512)


fig = plt.figure(figsize=(12,13))
grid = ImageGrid(fig, 111, nrows_ncols=(len(fields),4),
               axes_pad=0, label_mode='1', share_all=True,
               cbar_mode='edge', cbar_location='right',
               cbar_pad=0)

grid.axes_llc.tick_params(labelleft=False, labelbottom=False)
for ax in grid:
    ax.tick_params(which='both', axis='both', direction='in')
    ax.xaxis.set_major_locator(FixedLocator([-300,-200,-100,0,100,200,300]))

d_norm = LogNorm(1e-32, 1e-26)
t_norm = LogNorm(1e3, 1e8)
k_norm = LogNorm(1e0, 1e6)
p_norm = LogNorm(1e-18, 1e-14)
v_norm = SymLogNorm(1, linscale=0.2, base=10, vmin=-3e3, vmax=3e3)

ent_lo_cmap = plt.get_cmap('crest')(np.linspace(0, 1, 86))
ent_hi_cmap = plt.get_cmap('flare_r')(np.linspace(0, 1, 170))
colors = np.vstack((ent_lo_cmap, ent_hi_cmap))
ent_cmap = LinearSegmentedColormap.from_list(colors)

ax = grid.axes_column[0]
bar = AnchoredSizeBar(ax[0].transData, 100, "100 kpc", 3, 
                       label_top=True, color='white', frameon=False,
                       borderpad=1, size_vertical=5,
                       fontproperties={'size':'x-large',
                                       'weight':'bold'})
ax[0].add_artist(bar)

circle = Circle((0,0), 206, transform=ax[1].transData,
                edgecolor='white', fill=False, ls='--')
ax[1].add_artist(circle)

ax[0].text(0.04, 0.9, f"{ds_fid20.current_time.to('Gyr'):.2f}", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_fid20 = ax[0].imshow(np.array(frb_fid20['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_fid20 = ax[1].imshow(np.array(frb_fid20['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)
k_fid20 = ax[2].imshow(np.array(frb_fid20['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
p_fid20 = ax[3].imshow(np.array(frb_fid20['pressure']),
                         origin='lower', extent=extent,
                         cmap='cividis', norm=p_norm)
v_fid20 = ax[4].imshow(np.array(frb_fid20['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)

ax = grid.axes_column[1]
ax[0].text(0.04, 0.9, f"{ds_fid40.current_time.to('Gyr'):.2f}", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_fid40 = ax[0].imshow(np.array(frb_fid40['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_fid40 = ax[1].imshow(np.array(frb_fid40['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)
k_fid40 = ax[2].imshow(np.array(frb_fid40['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
p_fid40 = ax[3].imshow(np.array(frb_fid40['pressure']),
                         origin='lower', extent=extent,
                         cmap='cividis', norm=p_norm)
v_fid40 = ax[4].imshow(np.array(frb_fid40['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)

ax = grid.axes_column[2]
ax[0].text(0.04, 0.9, f"{ds_fid60.current_time.to('Gyr'):.2f}", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_fid60 = ax[0].imshow(np.array(frb_fid60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_fid60 = ax[1].imshow(np.array(frb_fid60['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)
k_fid60 = ax[2].imshow(np.array(frb_fid60['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
p_fid60 = ax[3].imshow(np.array(frb_fid60['pressure']),
                         origin='lower', extent=extent,
                         cmap='cividis', norm=p_norm)
v_fid60 = ax[4].imshow(np.array(frb_fid60['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)

ax = grid.axes_column[3]
ax[0].text(0.04, 0.9, f"{ds_fid80.current_time.to('Gyr'):.2f}", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

d_fid80 = ax[0].imshow(np.array(frb_fid80['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_fid80 = ax[1].imshow(np.array(frb_fid80['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)
k_fid80 = ax[2].imshow(np.array(frb_fid80['entropy']),
                         origin='lower', extent=extent,
                         cmap=ent_cmap, norm=k_norm)
p_fid80 = ax[3].imshow(np.array(frb_fid80['pressure']),
                         origin='lower', extent=extent,
                         cmap='cividis', norm=p_norm)
v_fid80 = ax[4].imshow(np.array(frb_fid80['radial_velocity'])/1e5, 
                         origin='lower', extent=extent,
                         cmap='coolwarm', norm=v_norm)

d_cb = fig.colorbar(d_fid80, cax=grid.cbar_axes[0], extend='both')
t_cb = fig.colorbar(t_fid80, cax=grid.cbar_axes[1], extend='both')
k_cb = fig.colorbar(k_fid80, cax=grid.cbar_axes[2], extend='both')
p_cb = fig.colorbar(p_fid80, cax=grid.cbar_axes[3], extend='both')
v_cb = fig.colorbar(v_fid80, cax=grid.cbar_axes[4], extend='both')

d_cb.set_ticks(FixedLocator([1e-31, 1e-29, 1e-27, 1e-25, 1e-32, 1e-30, 1e-28, 1e-26]))

t_cb.set_ticks(FixedLocator([1e3,1e4,1e5,1e6,1e7,1e8]))
t_cb.minorticks_off()

k_cb.set_ticks(FixedLocator([1e0, 1e2, 1e4, 1e6, 1e-1, 1e1, 1e3, 1e5]))

p_cb.set_ticks(FixedLocator([1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12]))

v_cb.set_ticks(FixedLocator([-1e3,-1e2,-1e1,0,1e1,1e2,1e3]))
v_cb.minorticks_off()

d_cb.set_label(r'Density  [g cm$^{-3}$]')
t_cb.set_label(r'Temperature  [K]')
p_cb.set_label(r'Pressure  [erg cm$^{-3}$]')
k_cb.set_label(r'Entropy  [keV cm$^2$]')
v_cb.set_label(r'Radial Velocity  [km/s]')

fig.subplots_adjust(left=0.02, right=0.95, bottom=0.01, top=0.99)
fig.savefig("fig_edge-ev_fid.pdf")

