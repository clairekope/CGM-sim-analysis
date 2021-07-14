
#!/usr/bin/env python
# coding: utf-8

#########################################################################
# Plot multipanel figures comparing sim variants (stored locally) @ 3 Gyr
# using face on view of disk
#########################################################################

import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator, MultipleLocator
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

ds_fid60 = yt.load("../sample_data/fid/DD0060/DD0060") # 3 Gyr
ds_cfw60 = yt.load("../sample_data/cflow/DD0060/DD0060")
ds_low60 = yt.load("../sample_data/tctff5/DD0060/DD0060")
ds_hih60 = yt.load("../sample_data/tctff20/DD0060/DD0060")
ds_lin60 = yt.load("../sample_data/linrot/DD0060/DD0060")
ds_nor60 = yt.load("../sample_data/norot/DD0060/DD0060")

# code length is the same in all sims
center = yt.YTQuantity(ds_fid60.quan(0.5,'code_length').to('cm'))
thickness = yt.YTQuantity(1.3,'kpc')
width = yt.YTQuantity(100, 'kpc')
extent = (-width/2, width/2, -width/2, width/2)

fields = ['density','temperature']

def prep_frb(ds):
    rect = ds.region([center, center, center],
                     [center-width/2, center-width/2, center-thickness/2],
                     [center+width/2, center+width/2, center+thickness/2])

    p = yt.SlicePlot(ds, 'z', fields, width=width,
                          data_source=rect)#, weight_field='ones')
                          
    frb = p.data_source.to_frb(width, 512)

    return frb

frb_fid60 = prep_frb(ds_fid60)
frb_cfw60 = prep_frb(ds_cfw60)
frb_low60 = prep_frb(ds_low60)
frb_hih60 = prep_frb(ds_hih60)
frb_lin60 = prep_frb(ds_lin60)
frb_nor60 = prep_frb(ds_nor60)

fig = plt.figure(figsize=(11,9))
grid = ImageGrid(fig, 111, nrows_ncols=(3,2*len(fields)),
               axes_pad=0, label_mode='1', share_all=True,
               cbar_mode='edge', cbar_location='top',
               cbar_pad=0)

grid.axes_llc.tick_params(labelleft=False, labelbottom=False)
for ax in grid:
    ax.tick_params(which='both', axis='both', direction='in')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    
d_norm = LogNorm(1e-30, 1e-23)
t_norm = LogNorm(1e3, 5e6)

ax = grid.axes_row[0]
bar = AnchoredSizeBar(ax[0].transData, 10, "10 kpc", 3, 
                       label_top=True, color='white', frameon=False,
                       borderpad=1, size_vertical=2,
                       fontproperties={'size':'x-large',
                                       'weight':'bold'})
ax[0].add_artist(bar)

circle = Circle((0,0), 6.874581203331226, transform=ax[0].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[0].add_artist(circle)
circle = Circle((0,0), 6.874581203331226, transform=ax[1].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[1].add_artist(circle)

ax[0].text(0.04, 0.9, "Fiducial", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large',
                     'weight':'bold',
                     'color':'white'})

d_fid60 = ax[0].imshow(np.array(frb_fid60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_fid60 = ax[1].imshow(np.array(frb_fid60['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)

ax = grid.axes_row[1]
ax[0].text(0.04, 0.9, "CoolFlow", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

circle = Circle((0,0), 9.607565846320332, transform=ax[0].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[0].add_artist(circle)
circle = Circle((0,0), 9.607565846320332, transform=ax[1].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[1].add_artist(circle)

d_cfw60 = ax[0].imshow(np.array(frb_cfw60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_cfw60 = ax[1].imshow(np.array(frb_cfw60['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)

ax = grid.axes_row[2]
ax[0].text(0.04, 0.9, "LowRatio", 
           transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

circle = Circle((0,0), 10.312465486062305, transform=ax[0].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[0].add_artist(circle)
circle = Circle((0,0), 10.312465486062305, transform=ax[1].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[1].add_artist(circle)

d_low60 = ax[0].imshow(np.array(frb_low60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_low60 = ax[1].imshow(np.array(frb_low60['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)

ax = grid.axes_row[0]
ax[2].text(0.04, 0.9, "HighRatio", 
           transform=ax[2].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

circle = Circle((0,0), 5.445915564113368, transform=ax[2].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[2].add_artist(circle)
circle = Circle((0,0), 5.445915564113368, transform=ax[3].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[3].add_artist(circle)

d_hih60 = ax[2].imshow(np.array(frb_hih60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_hih60 = ax[3].imshow(np.array(frb_hih60['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)
                         
ax = grid.axes_row[1]
ax[2].text(0.04, 0.9, "LinRot", 
           transform=ax[2].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

circle = Circle((0,0), 8.578234112925314, transform=ax[2].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[2].add_artist(circle)
circle = Circle((0,0), 8.578234112925314, transform=ax[3].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[3].add_artist(circle)

d_lin60 = ax[2].imshow(np.array(frb_lin60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_lin60 = ax[3].imshow(np.array(frb_lin60['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)

ax = grid.axes_row[2]
ax[2].text(0.04, 0.9, "NoRot", 
           transform=ax[2].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

circle = Circle((0,0), 7.236936862951981, transform=ax[2].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[2].add_artist(circle)
circle = Circle((0,0), 7.236936862951981, transform=ax[3].transData,
                edgecolor='white', fill=False, ls='-', lw=1)
ax[3].add_artist(circle)

d_nor60 = ax[2].imshow(np.array(frb_nor60['density']),
                         origin='lower', extent=extent,
                         norm=d_norm)
t_nor60 = ax[3].imshow(np.array(frb_nor60['temperature']),
                         origin='lower', extent=extent,
                         cmap='magma', norm=t_norm)

d_cb = fig.colorbar(d_fid60, cax=grid.cbar_axes[0], extend='both', orientation='horizontal', ticklocation='top')
t_cb = fig.colorbar(t_fid60, cax=grid.cbar_axes[1], extend='both', orientation='horizontal', ticklocation='top')

d_cb2 = fig.colorbar(d_fid60, cax=grid.cbar_axes[2], extend='both', orientation='horizontal', ticklocation='top')
t_cb2 = fig.colorbar(t_fid60, cax=grid.cbar_axes[3], extend='both', orientation='horizontal', ticklocation='top')

d_cb.ax.xaxis.set_minor_locator(FixedLocator([1e-29, 1e-27, 1e-25, 1e-23]))
d_cb.ax.xaxis.set_major_locator(FixedLocator([1e-30, 1e-28, 1e-26, 1e-24]))
d_cb.ax.xaxis.set_minor_formatter(NullFormatter())
d_cb2.ax.xaxis.set_minor_locator(FixedLocator([1e-29, 1e-27, 1e-25, 1e-23]))
d_cb2.ax.xaxis.set_major_locator(FixedLocator([1e-30, 1e-28, 1e-26, 1e-24]))
d_cb2.ax.xaxis.set_minor_formatter(NullFormatter())

t_cb.set_ticks(FixedLocator([1e3,1e4,1e5,1e6,1e7,1e8]))
#t_cb.minorticks_off()
t_cb2.set_ticks(FixedLocator([1e3,1e4,1e5,1e6,1e7,1e8]))
#t_cb2.minorticks_off()

d_cb.set_label(r'Density  [g cm$^{-3}$]', fontsize='large')
t_cb.set_label(r'Temperature  [K]', fontsize='large')
d_cb2.set_label(r'Density  [g cm$^{-3}$]', fontsize='large')
t_cb2.set_label(r'Temperature  [K]', fontsize='large')

fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.95)
fig.savefig("../fig_face-comp_slices.eps")

