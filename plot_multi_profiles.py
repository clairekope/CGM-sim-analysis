#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_profiles_over_time(xind, yind, 
                            tstart, tend, dt,
                            dtDD=5, cmap='cividis', cmap_label='t (Myr)',
                            fpath='./', fbase='profiles', fpad=4, fext='txt'):
    """
    Plot profiles designated by columns <xind> and <yind>
    in files <fpath>/<fbase>_{0<fpad>d}.<fext>
    at times np.arange(<tstart>, <tend>, dt)
    using the colormap <cmap>.
    Times are in code units and converted to Enzo datadumps via <dtDD>
    """
                            
    times = np.arange(tstart, tend+dt, dt, dtype=int)
    ntimes = times.size

    c = np.linspace(tstart, tend, ntimes)
    cmap = mpl.cm.get_cmap(cmap, ntimes)
    norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
    cm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cm.set_array([])

    fig, ax = plt.subplots(dpi=400)

    for n in times//dtDD:
        
        f = fpath + fbase + '_' + str(n).zfill(fpad) + '.' + fext
        x, y = np.genfromtxt(f, usecols=(xind, yind), unpack=True)
    
        ax.loglog(x, y, c=cm.to_rgba((n-1)*10))

    fig.colorbar(cm, label=cmap_label)
                 #ticks=mpl.ticker.LinearLocator(times[0]-dt/2, times[-1]+dt/2))
    
    return fig, ax

