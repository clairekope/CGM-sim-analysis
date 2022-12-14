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

        fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(23,15))

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


        #
        # Phase Diagrams
        #
        dt_ph = yt.PhasePlot(rect, 'density', 'temperature', 'cell_mass').profile
        pk_ph = yt.PhasePlot(rect, 'pressure', 'entropy', 'cell_mass').profile
        rk_ph = yt.PhasePlot(rect, 'radius', 'entropy', 'cell_mass').profile

        # 16th, 50th, and 84th precentile 1D profiles
        dt_med = np.ones(dt_ph.x_bins.size-1) * np.nan
        dt_16 = np.ones(dt_ph.x_bins.size-1) * np.nan
        dt_84 = np.ones(dt_ph.x_bins.size-1) * np.nan
        
        d_binner = np.digitize(rect['density'], dt_ph.x_bins)
        for i in range(1, dt_ph.x_bins.size):
            this_bin = d_binner==i
            try:
                dt_med[i-1] = wq.median(rect['temperature'][this_bin],
                                        rect['cell_mass'][this_bin])
                dt_16[i-1] = wq.quantile(rect['temperature'][this_bin],
                                         rect['cell_mass'][this_bin],
                                         0.16)
                dt_84[i-1] = wq.quantile(rect['temperature'][this_bin],
                                         rect['cell_mass'][this_bin],
                                         0.84)
            except ValueError:
                continue
        
        pk_med = np.ones(pk_ph.x_bins.size-1) * np.nan
        pk_16 = np.ones(pk_ph.x_bins.size-1) * np.nan
        pk_84 = np.ones(pk_ph.x_bins.size-1) * np.nan
        
        p_binner = np.digitize(rect['pressure'], pk_ph.x_bins)
        for i in range(1, pk_ph.x_bins.size):
            this_bin = p_binner==i
            try:
                pk_med[i-1] = wq.median(rect['entropy'][this_bin],
                                        rect['cell_mass'][this_bin])
                pk_16[i-1] = wq.quantile(rect['entropy'][this_bin],
                                         rect['cell_mass'][this_bin],
                                         0.16)
                pk_84[i-1] = wq.quantile(rect['entropy'][this_bin],
                                         rect['cell_mass'][this_bin],
                                         0.84)
            except ValueError:
                continue

        rk_med = np.ones(rk_ph.x_bins.size-1) * np.nan
        rk_16 = np.ones(rk_ph.x_bins.size-1) * np.nan
        rk_84 = np.ones(rk_ph.x_bins.size-1) * np.nan
        
        r_binner = np.digitize(rect['radius'], rk_ph.x_bins)
        for i in range(1, rk_ph.x_bins.size):
            this_bin = r_binner==i
            try:
                rk_med[i-1] = wq.median(rect['entropy'][this_bin],
                                        rect['cell_mass'][this_bin])
                rk_16[i-1] = wq.quantile(rect['entropy'][this_bin],
                                         rect['cell_mass'][this_bin],
                                         0.16)
                rk_84[i-1] = wq.quantile(rect['entropy'][this_bin],
                                         rect['cell_mass'][this_bin],
                                         0.84)
            except ValueError:
                continue
        
        # Plot phase diagrams
        dt_im = ax[2,0].pcolormesh(dt_ph.x_bins, dt_ph.y_bins,
                                   dt_ph['cell_mass'].T.to('Msun'),
                                   norm=LogNorm(1e-4,1e6), cmap='viridis')
        
        pk_im = ax[2,1].pcolormesh(pk_ph.x_bins, pk_ph.y_bins,
                                   pk_ph['cell_mass'].T.to('Msun'),
                                   norm=LogNorm(1e-4,1e6), cmap='cividis')

        rk_im = ax[2,2].pcolormesh(rk_ph.x_bins.to('kpc'), rk_ph.y_bins,
                                   rk_ph['cell_mass'].T.to('Msun'),
                                   norm=LogNorm(1e-3,1e7), cmap='magma')

        ax[2,0].plot(dt_ph.x, dt_med, 'k-')
        ax[2,0].plot(dt_ph.x, dt_16, 'k--')
        ax[2,0].plot(dt_ph.x, dt_84, 'k--')

        ax[2,1].plot(pk_ph.x, pk_med, 'k-')
        ax[2,1].plot(pk_ph.x, pk_16, 'k--')
        ax[2,1].plot(pk_ph.x, pk_84, 'k--')

        ax[2,2].plot(rk_ph.x.to('kpc'), rk_med, 'k-')
        ax[2,2].plot(rk_ph.x.to('kpc'), rk_16, 'k--')
        ax[2,2].plot(rk_ph.x.to('kpc'), rk_84, 'k--')

        # Adjust limits
        ax[2,0].set_xlim(1e-34, 1e-21)
        ax[2,0].set_xlabel(r'Density [g cm$^{-3}$]')
        ax[2,0].set_ylim(1e1, 1e9)
        ax[2,0].set_ylabel('Temperature [K]')

        ax[2,1].set_xlim(1e-22, 1e-10)
        ax[2,1].set_xlabel(r'Pressure [dyn cm$^{-2}$]')
        ax[2,1].set_ylim(1e-7, 1e7)
        ax[2,1].set_ylabel(r'Entropy [keV cm$^2$]')

        ax[2,2].set_xlim(1e-1, 1e3)
        ax[2,2].set_xlabel('Radius [kpc]')
        ax[2,2].set_ylim(1e-7, 1e7)
        ax[2,2].set_ylabel(r'Entropy [keV cm$^2$]')

        # Add colorbars and labels
        dt_cb = fig.colorbar(dt_im, ax=ax[2,0], pad=0.01, shrink=0.9)
        pk_cb = fig.colorbar(pk_im, ax=ax[2,1], pad=0.01, shrink=0.9)
        rk_cb = fig.colorbar(rk_im, ax=ax[2,2], pad=0.01, shrink=0.9)

        dt_cb.set_label(r'Cell Mass [M$_\odot$]')
        pk_cb.set_label(r'Cell Mass [M$_\odot$]')
        rk_cb.set_label(r'Cell Mass [M$_\odot$]')

        # Set plot view stuff
        for i in range(3):
            ax[2,i].set_xscale('log')
            ax[2,i].set_yscale('log')
            ax[2,i].set_aspect( 1.0 / ax[2,i].get_data_ratio() )


        #
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
        fig.savefig('{}_thermo_slab_proj_min_{}.png'.format(ds.basename, label))

        plt.close(fig)
