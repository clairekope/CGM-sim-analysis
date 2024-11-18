#!/usr/bin/env python
# coding: utf-8

import yt
yt.enable_parallelism()
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FixedLocator, MultipleLocator
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


datasets = yt.load('DD????/DD????')
for ds in datasets.piter(dynamic=False, ):

    center = ds.quan(0.5,'code_length').to('cm')
    thickness = ds.quan(3.5,'kpc')
    width_edge = ds.quan(400, 'kpc')
    width_face = ds.quan(100, 'kpc')

    extent_edge = (-width_edge/2, width_edge/2, -width_edge/2, width_edge/2)
    #extent_face = (-width_face/2, width_face/2, -width_face/2, width_face/2)

    rect_edge = ds.region([center, center, center],
                     [center-thickness/2, center-width_edge/2, center-width_edge/2],
                     [center+thickness/2, center+width_edge/2, center+width_edge/2])

    rect_face = ds.region([center, center, center],
                     [center-width_face/2, center-width_face/2, center-thickness/2],
                     [center+width_face/2, center+width_face/2, center+thickness/2])

    p_edge = yt.ProjectionPlot(ds, 'x', ['density','temperature'], weight_field='ones',
                               width=width_edge, data_source=rect_edge)


    p_face = yt.ProjectionPlot(ds, 'z', ['density','temperature'], weight_field='ones',
                               width=width_face, data_source=rect_face)


    frb_edge = p_edge.data_source.to_frb(width_edge, 512)
    frb_face = p_face.data_source.to_frb(width_face, 512)


    fig = plt.figure(figsize=(7,6))
    grid = ImageGrid(fig, 111, nrows_ncols=(2,2), aspect=True,
                   axes_pad=0, label_mode='1', share_all=True,
                   cbar_mode='edge', cbar_location='right', cbar_pad=0)

    grid.axes_llc.tick_params(labelleft=False, labelbottom=False)
    for ax in grid:
        ax.tick_params(which='both', axis='both', direction='in')
        
    grid[0].xaxis.set_major_locator(MultipleLocator(40))
    grid[0].xaxis.set_minor_locator(MultipleLocator(20))
    grid[2].xaxis.set_major_locator(MultipleLocator(40))
    grid[2].xaxis.set_minor_locator(MultipleLocator(20))

    grid[1].yaxis.set_major_locator(MultipleLocator(50))
    grid[1].yaxis.set_minor_locator(MultipleLocator(25))
    grid[3].yaxis.set_major_locator(MultipleLocator(50))
    grid[3].yaxis.set_minor_locator(MultipleLocator(25))

    d_norm = LogNorm(1e-32, 1e-24)
    t_norm = LogNorm(1e3, 1e8)

    ax = grid.axes_row[0]
    d_face = ax[0].imshow(np.array(frb_face['density']),
                          origin='lower', extent=extent_edge,
                          norm=d_norm)
    d_edge = ax[1].imshow(np.array(frb_edge['density']),
                          origin='lower', extent=extent_edge,
                          norm=d_norm)
    
    ax[0].text(0.04, 0.9, f"{ds.current_time.to('Gyr'):.2f}", transform=ax[0].transAxes,
           fontdict={'size':'x-large','weight':'bold','color':'white'})

    ax = grid.axes_row[1]
    t_face = ax[0].imshow(np.array(frb_face['temperature']),
                          origin='lower', extent=extent_edge,
                          cmap='magma', norm=t_norm)
    t_edge = ax[1].imshow(np.array(frb_edge['temperature']),
                          origin='lower', extent=extent_edge,
                          cmap='magma', norm=t_norm)

    # face's pixel scale is inflated by a factor of 4
    # so 40 kpc of axis size is 10 kpc of physical size
    bar_face = AnchoredSizeBar(ax[1].transData, 40, "10 kpc", 3, 
                               label_top=True, color='white', frameon=False,
                               borderpad=1, size_vertical=5,
                               fontproperties={'size':'x-large',
                                               'weight':'bold'})
    bar_edge = AnchoredSizeBar(ax[1].transData, 100, "100 kpc", 3, 
                               label_top=True, color='white', frameon=False,
                               borderpad=1, size_vertical=5,
                               fontproperties={'size':'x-large',
                                               'weight':'bold'})
    ax[0].add_artist(bar_face)
    ax[1].add_artist(bar_edge)

    d_cb = fig.colorbar(d_edge, cax=grid.cbar_axes[0], extend='both')
    t_cb = fig.colorbar(t_edge, cax=grid.cbar_axes[1], extend='both')

    d_cb.set_ticks(FixedLocator([1e-32,1e-31,1e-30,1e-29,1e-28,1e-27,1e-26,1e-25,1e-24]))

    d_cb.set_label(r'Density  [g cm$^{-3}$]', fontsize='x-large')
    t_cb.set_label(r'Temperature  [K]', fontsize='x-large')

    fig.subplots_adjust(left=0.02, right=0.87, bottom=0.01, top=0.99)
    fig.savefig(f"{ds.basename}_4panel.png", dpi=300)
