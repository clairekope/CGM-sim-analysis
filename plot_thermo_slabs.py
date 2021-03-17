import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
from scipy.stats import binned_statistic
from calc_enclosed_mass import *

rcParams.update({'font.size': 14})

masses = np.genfromtxt('masses_over_time.txt')
sfr = np.genfromtxt('sfr.txt')

datasets = yt.load('DD????/DD????')
for ds in datasets.piter(dynamic=False, ):

    center = ds.quan(0.5,'code_length')
    rs = ds.quan(3.5,'kpc')
    
    widths = [
        ds.quan(100,'kpc'),
        ds.quan(100,'kpc'),
        ds.quan(200,'kpc'), 
        ds.quan(400,'kpc'),
        # ds.quan(800,'kpc')
        ]
    thicknesses = [
        rs,
        rs,
        rs, 
        rs,
        # rs*2
        ]
    labels = [
        '100kpc',
        '100kpc_face',
        '200kpc',
        '400kpc',
        # '800kpc'
        ]

    for width, thickness, label in zip(widths, thicknesses, labels):

        fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(24.5,15.5))

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

        # Sphere is used for calculating average tff
        sph = ds.sphere([center, center, center], width)

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
    
        d_im = ax[0,0].imshow(d_arr, origin='lower', norm=LogNorm(1e-32,1e-24), 
                              extent=extent)

        t_im = ax[1,0].imshow(t_arr, origin='lower', norm=LogNorm(1e3,1e8),
                              extent=extent, cmap='magma')

        p_im = ax[0,1].imshow(p_arr, origin='lower', norm=LogNorm(1e-21,1e-13),
                              extent=extent, cmap='inferno')

        k_im = ax[1,1].imshow(k_arr, origin='lower', norm=LogNorm(1e-2,1e6),
                              extent=extent, cmap='cividis')

        v_im = ax[0,2].imshow(v_arr, origin='lower',
                              norm=SymLogNorm(1, linscale=0.3, base=10,
                                              vmin=-3e3, vmax=3e3),
                              extent=extent, cmap='coolwarm')

        c_im = ax[1,2].imshow(c_arr, origin='lower', norm=LogNorm(1e-2,1e4),
                              extent=extent, cmap='twilight')

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
        rc_ph = yt.PhasePlot(rect, 'radius', 'cooling_time', 'cell_mass').profile
        rk_ph = yt.PhasePlot(rect, 'radius', 'entropy', 'cell_mass').profile

        # 16th, 50th, and 84th precentile 1D profiles
        # e.g., median temperature in each density bin,
        # or median entropy in each radial bin
        dt_med = binned_statistic(rect['density'], rect['temperature'],
                                  statistic='median', bins=dt_ph.x_bins)
        dt_16 = binned_statistic(rect['density'], rect['temperature'],
                                 statistic=lambda y: np.percentile(y,16),
                                 bins=dt_ph.x_bins)
        dt_84 = binned_statistic(rect['density'], rect['temperature'],
                                 statistic=lambda y: np.percentile(y,84),
                                 bins=dt_ph.x_bins)

        pk_med = binned_statistic(rect['pressure'], rect['entropy'],
                                  statistic='median', bins=pk_ph.x_bins)
        pk_16 = binned_statistic(rect['pressure'], rect['entropy'],
                                 statistic=lambda y: np.percentile(y,16),
                                 bins=pk_ph.x_bins)
        pk_84 = binned_statistic(rect['pressure'], rect['entropy'],
                                 statistic=lambda y: np.percentile(y,84),
                                 bins=pk_ph.x_bins)

        rc_med = binned_statistic(rect['radius'].to('kpc'),
                                  rect['cooling_time'].to('Gyr'),
                                  statistic='median',
                                  bins=rc_ph.x_bins.to('kpc'))
        rc_16 = binned_statistic(rect['radius'].to('kpc'),
                                 rect['cooling_time'].to('Gyr'),
                                 statistic=lambda y: np.percentile(y,16),
                                 bins=rc_ph.x_bins.to('kpc'))
        rc_84 = binned_statistic(rect['radius'].to('kpc'),
                                 rect['cooling_time'].to('Gyr'),
                                 statistic=lambda y: np.percentile(y,84),
                                 bins=rc_ph.x_bins.to('kpc'))

        rk_med = binned_statistic(rect['radius'].to('kpc'), rect['entropy'],
                                  statistic='median', bins=rk_ph.x_bins.to('kpc'))
        rk_16 = binned_statistic(rect['radius'].to('kpc'), rect['entropy'],
                                 statistic=lambda y: np.percentile(y,16),
                                 bins=rk_ph.x_bins.to('kpc'))
        rk_84 = binned_statistic(rect['radius'].to('kpc'), rect['entropy'],
                                 statistic=lambda y: np.percentile(y,84),
                                 bins=rk_ph.x_bins.to('kpc'))

        # Calculate tff for tcool vs tff plot
        # tff is calculated for material enclosed within a sphere
        # as MN_accel averages the stellar potential over samples of a sphere
        # (technically of a circle b/c of azimuthal symmetry)
        Mtot_r = NFW_mass_enclosed(rc_ph.x) + \
                 cell_and_particle_mass_enclosed(sph, rc_ph.x_bins.v)

        g_r = yt.units.G * Mtot_r / np.power(rc_ph.x, 2)
        g_r = g_r.to('cm/s**2') + MN_accel(rc_ph.x)

        t_ff = np.sqrt( 2*rc_ph.x / g_r )
        
        # Plot phase diagrams
        dt_im = ax[2,0].pcolormesh(dt_ph.x_bins, dt_ph.y_bins,
                                   dt_ph['cell_mass'].T.to('Msun'),
                                   norm=LogNorm(1e-4,1e6), cmap='viridis')
        
        pk_im = ax[2,1].pcolormesh(pk_ph.x_bins, pk_ph.y_bins,
                                   pk_ph['cell_mass'].T.to('Msun'),
                                   norm=LogNorm(1e-4,1e6), cmap='cividis')

        rc_im = ax[2,2].pcolormesh(rc_ph.x_bins.to('kpc'), rc_ph.y_bins.to('Gyr'),
                                   rc_ph['cell_mass'].T.to('Msun'),
                                   norm=LogNorm(1e-3,1e7), cmap='inferno')

        rk_im = ax[2,3].pcolormesh(rk_ph.x_bins.to('kpc'), rk_ph.y_bins,
                                   rk_ph['cell_mass'].T.to('Msun'),
                                   norm=LogNorm(1e-3,1e7), cmap='magma')

        # Plot percentiles
        ax[2,0].plot(dt_ph.x, dt_med[0], 'k-')
        ax[2,0].plot(dt_ph.x, dt_16[0], 'k--')
        ax[2,0].plot(dt_ph.x, dt_84[0], 'k--')

        ax[2,1].plot(pk_ph.x, pk_med[0], 'k-')
        ax[2,1].plot(pk_ph.x, pk_16[0], 'k--')
        ax[2,1].plot(pk_ph.x, pk_84[0], 'k--')

        ax[2,2].plot(rc_ph.x.to('kpc'), rc_med[0], 'w-')
        ax[2,2].plot(rc_ph.x.to('kpc'), rc_16[0], 'w--')
        ax[2,2].plot(rc_ph.x.to('kpc'), rc_84[0], 'w--')

        ax[2,3].plot(rk_ph.x.to('kpc'), rk_med[0], 'w-')
        ax[2,3].plot(rk_ph.x.to('kpc'), rk_16[0], 'w--')
        ax[2,3].plot(rk_ph.x.to('kpc'), rk_84[0], 'w--')

        # Plot tff
        ax[2,2].plot(rc_ph.x.to('kpc'), t_ff.to('Gyr'), 'g-')
        ax[2,2].plot(rc_ph.x.to('kpc'), 10*t_ff.to('Gyr'), 'g-.')

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
        ax[2,2].set_ylim(1e-6, 1e8)
        ax[2,2].set_ylabel('Cooling Time [Gyr]')

        ax[2,3].set_xlim(1e-1, 1e3)
        ax[2,3].set_xlabel('Radius [kpc]')
        ax[2,3].set_ylim(1e-7, 1e7)
        ax[2,3].set_ylabel(r'Entropy [keV cm$^2$]')

        # Add colorbars and labels
        dt_cb = fig.colorbar(dt_im, ax=ax[2,0], pad=0.01, shrink=0.9)
        pk_cb = fig.colorbar(pk_im, ax=ax[2,1], pad=0.01, shrink=0.9)
        rc_cb = fig.colorbar(rc_im, ax=ax[2,2], pad=0.01, shrink=0.9)
        rk_cb = fig.colorbar(rk_im, ax=ax[2,3], pad=0.01, shrink=0.9)

        dt_cb.set_label(r'Cell Mass [M$_\odot$]')
        pk_cb.set_label(r'Cell Mass [M$_\odot$]')
        rc_cb.set_label(r'Cell Mass [M$_\odot$]')
        rk_cb.set_label(r'Cell Mass [M$_\odot$]')

        # Set plot view stuff
        for i in range(4):
            ax[2,i].set_xscale('log')
            ax[2,i].set_yscale('log')
            ax[2,i].set_aspect( 1.0 / ax[2,i].get_data_ratio() )

        #
        # Plots over time
        #

        # Translate output number to index
        index = int(ds.basename[-4:])

        # Plot SFR up to this point
        time = masses[index,0]
        past_sfr = sfr[:,0] <= time
        ax[0,3].plot(sfr[past_sfr,0], sfr[past_sfr, 1], 'k')
        ax[0,3].set_xlim(0, 4000)
        ax[0,3].set_xlabel('Time [Myr]')
        ax[0,3].set_ylim(0, 47)
        ax[0,3].set_ylabel(r'SFR [M$_\odot$/yr]')

        # Plot sum of SFR up to this point
        sum_ax = ax[0,3].twinx()
        sum_ax.plot(sfr[past_sfr,0], np.nancumsum(sfr[past_sfr, 1]),
                    color='C0', ls='--')
        sum_ax.set_ylim(0, 10000)
        sum_ax.tick_params(axis='y', labelcolor='C0')
        sum_ax.set_ylabel(r'$\sum$SFR [M$_\odot$/yr]', color='C0')
        
        ax[0,3].xaxis.tick_top()
        ax[0,3].xaxis.set_label_position('top')
        ax[0,3].yaxis.set_major_locator(MultipleLocator(5))

        # Plot disk masses up to this point
        ax[1,3].plot(masses[:index,0], masses[:index,1],
                     ls='--', label=r'$\rm M_{disk}$')
        ax[1,3].plot(masses[:index,0], masses[:index,2],
                     ls='-.', label=r'$\rm M_\ast$')
        ax[1,3].plot(masses[:index,0], masses[:index,3],
                     ls='-', label=r'$\rm M_{sum}$')       
        ax[1,3].set_xlim(0, 4000)
        ax[1,3].set_xlabel('Time [Myr]')
        ax[1,3].set_ylim(0, 1e10)
        ax[1,3].set_ylabel(r'Mass [M$_\odot$]')
        ax[1,3].legend()

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
        fig.savefig('{}_thermo_slab_proj_{}.png'.format(ds.basename, label))

        plt.close(fig)
