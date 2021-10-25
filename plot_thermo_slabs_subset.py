####################################################
# Plot a large 12 panel plot with thin slab
# projection, phase plots, SFR, and disk mass growth
####################################################

import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np

rcParams.update({'font.size': 14})

datasets = yt.load('DD00[0-2]?/DD00??')
for ds in datasets.piter(dynamic=False, ):

    center = ds.quan(0.5,'code_length')
    rs = ds.quan(3.5,'kpc')
    
    widths = [
        ds.quan(100,'kpc'),
        # ds.quan(100,'kpc'),
        ds.quan(200,'kpc'), 
        ds.quan(400,'kpc'),
        ds.quan(800,'kpc')
        ]
    thicknesses = [
        rs,
        # rs,
        rs, 
        rs,
        rs*2
        ]
    labels = [
        '100kpc',
        # '100kpc_face',
        '200kpc',
        '400kpc',
        '800kpc'
        ]

    for width, thickness, label in zip(widths, thicknesses, labels):

        fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(22,15))

        if 'face' in label:
            view = 'z'

            # Thin slab is used for projections so I can avoid weighting by e.g. mass
            rect = ds.region([center, center, center],
                             [center-width/2, center-width/2, center-thickness/2],
                             [center+width/2, center+width/2, center+thickness/2])
            
        else:
            view = 'x'
            
            # Thin slab is used for projections so I can avoid weighting by e.g. mass
            rect = ds.region([center, center, center],
                             [center-thickness/2, center-width/2, center-width/2],
                             [center+thickness/2, center+width/2, center+width/2])

        assert (rect.get_field_parameter("center") == center).all()

        #
        # Projection Plots
        #
        fields = ['density','temperature','pressure','entropy',
                  'radial_velocity','cooling_time']

        p = yt.ProjectionPlot(ds, view, fields, width=width,
                              data_source=rect, weight_field ='ones')

        # Extract fixed resolution buffers
        p_frb = p.data_source.to_frb(width, 512)

        d_arr = np.array(p_frb['density'])
        t_arr = np.array(p_frb['temperature'])
        p_arr = np.array(p_frb['pressure'])
        k_arr = np.array(p_frb['entropy'])
        v_arr = np.array(p_frb['radial_velocity']) / 1e5 # km/s
        c_arr = np.array(p_frb['cooling_time']) / 3.16887e8 / 1e9 # Gyr

        # Plot FRBs
        extent = (-width/2, width/2, -width/2, width/2)
    
        d_im = ax[0,0].imshow(d_arr, origin='lower', norm=LogNorm(1e-31,1e-23), 
                              extent=extent)

        t_im = ax[1,0].imshow(t_arr, origin='lower', norm=LogNorm(1e3,1e7),
                              extent=extent, cmap='magma')

        p_im = ax[0,1].imshow(p_arr, origin='lower', norm=LogNorm(1e-20,1e-12),
                              extent=extent, cmap='inferno')

        k_im = ax[1,1].imshow(k_arr, origin='lower', norm=LogNorm(1e-2,1e3),
                              extent=extent, cmap='cividis')

        v_im = ax[0,2].imshow(v_arr, origin='lower',
                              norm=SymLogNorm(1, linscale=0.3, base=10,
                                              vmin=-2e2, vmax=2e2),
                              extent=extent, cmap='coolwarm')

        c_im = ax[1,2].imshow(c_arr, origin='lower', norm=LogNorm(1e-3,1e2),
                              extent=extent, cmap='plasma')

        # Add colorbars
        d_cb = fig.colorbar(d_im, ax=ax[0,0], pad=0.01, shrink=0.9)
        t_cb = fig.colorbar(t_im, ax=ax[1,0], pad=0.01, shrink=0.9)
        p_cb = fig.colorbar(p_im, ax=ax[0,1], pad=0.01, shrink=0.9)
        k_cb = fig.colorbar(k_im, ax=ax[1,1], pad=0.01, shrink=0.9)
        v_cb = fig.colorbar(v_im, ax=ax[0,2], pad=0.01, shrink=0.9)
        c_cb = fig.colorbar(c_im, ax=ax[1,2], pad=0.01, shrink=0.9)

        # Add colorbar labels
        d_cb.set_label(r'< Density > [g cm$^{-3}$]')
        t_cb.set_label(r'< Temperature > [K]')
        p_cb.set_label(r'< Pressure > [dyn cm$^{-2}$]')
        k_cb.set_label(r'< Entropy > [keV cm$^2$]')
        v_cb.set_label(r'< Radial Velocity > [km/s]')
        c_cb.set_label(r'< Cooling Time > [Gyr]')

       
        # Final details & save
        #
        if 'face' in label:
            ax[0,0].set_xlabel('x (kpc)')
            ax[0,0].set_ylabel('y (kpc)')
        else:
            ax[0,0].set_xlabel('y (kpc)')
            ax[0,0].set_ylabel('z (kpc)')
        
        time = '{:.0f} Myr'.format(np.round(ds.current_time.to('Myr'),0))
        ax[0,0].text(0.03, 0.03, time, transform=ax[0,0].transAxes,
                     color='white', size='x-large', weight='bold')

        fig.suptitle("{} width, {} thickness".format(width, thickness), fontsize=28)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        fig.subplots_adjust(hspace=0.1, wspace=0.3)
        fig.savefig('{}_thermo_slab_proj_{}.png'.format(ds.basename, label))

        plt.close(fig)
