import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

datasets = yt.load('DD????/DD????')
for ds in datasets.piter():

    fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,
                           figsize=(12,10))

    center = ds.quan(0.5,'code_length')
    rs = ds.quan(3.5,'kpc')
    width = ds.quan(600,'kpc')
    rect = ds.region([center, center, center],
                     [center-rs/2, center-width/2, center-width/2],
                     [center+rs/2, center+width/2, center+width/2])

    fields = ['density','temperature','pressure','entropy']
    p = yt.ProjectionPlot(ds, 'x', fields,
                          width=width, data_source=rect, weight_field ='ones')

    p_frb = p.data_source.to_frb(width, 512)

    d_arr = np.array(p_frb['density'])
    t_arr = np.array(p_frb['temperature'])
    p_arr = np.array(p_frb['pressure'])
    k_arr = np.array(p_frb['entropy'])

    extent = (-width/2, width/2, -width/2, width/2)
    
    d_im = ax[0,0].imshow(d_arr, origin='lower', norm=LogNorm(1e-31,1e-24), 
                          extent=extent)

    t_im = ax[1,0].imshow(t_arr, origin='lower', norm=LogNorm(1e3,1e8),
                          extent=extent, cmap='magma')

    p_im = ax[0,1].imshow(p_arr, origin='lower', norm=LogNorm(1e-21,1e-14),
                          extent=extent, cmap='inferno')

    k_im = ax[1,1].imshow(k_arr, origin='lower', norm=LogNorm(1e-2,1e6),
                          extent=extent, cmap='cividis')

    d_cb = fig.colorbar(d_im, ax=ax[0,0], pad=0.01)
    t_cb = fig.colorbar(t_im, ax=ax[1,0], pad=0.01)
    p_cb = fig.colorbar(p_im, ax=ax[0,1], pad=0.01)
    k_cb = fig.colorbar(k_im, ax=ax[1,1], pad=0.01)

    d_cb.set_label(r'Density (g cm$^{-3}$)')
    t_cb.set_label(r'Temperature (K)')
    p_cb.set_label(r'Pressure (dyn cm$^{-2}$)')
    k_cb.set_label(r'Entropy (keV cm$^2$)')

    ax[1,0].set_xlabel('y (kpc)')
    ax[1,0].set_ylabel('z (kpc)')

    fig.tight_layout()
    fig.savefig('{}_thermo_slab_proj.png'.format(ds.basename))
