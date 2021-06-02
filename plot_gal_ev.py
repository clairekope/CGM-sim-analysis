#####################################################
# Save multipanel plot with
# thin slab projections of density, temperature,
# and radial velocity with the disk edge on & face on
#####################################################

import yt
yt.enable_parallelism()
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import MultipleLocator
import numpy as np

rcParams.update({'font.size': 14})

datasets = yt.load('DD????/DD????')
for ds in datasets.piter(dynamic=False, ):
    center = ds.quan(0.5,'code_length')
    thickness = ds.quan(3.5,'kpc') # 1 rs
    width = ds.quan(100,'kpc')
    
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15,11),
                           gridspec_kw={"height_ratios":(0.05,1,1)})
    
    rect_edge = ds.region([center, center, center],
                          [center-width/2, center-thickness/2, center-width/2],
                          [center+width/2, center+thickness/2, center+width/2])
    rect_face = ds.region([center, center, center],
                          [center-width/2, center-width/2, center-thickness/2],
                          [center+width/2, center+width/2, center+thickness/2])
    
    fields = ['density', 'temperature', 'radial_velocity']
    
    p_edge = yt.ProjectionPlot(ds, 'y', fields, width=width,
                               data_source=rect_edge, weight_field="ones"
                              ).data_source.to_frb(width, 512)
    
    p_face = yt.ProjectionPlot(ds, 'z', fields, width=width,
                               data_source=rect_face, weight_field="ones"
                              ).data_source.to_frb(width, 512)
    
    d_edge = np.array(p_edge["density"]).T
    t_edge = np.array(p_edge["temperature"]).T
    v_edge = np.array(p_edge["radial_velocity"]).T / 1e5 # km/s
    d_face = np.array(p_face["density"])
    t_face = np.array(p_face["temperature"])
    v_face = np.array(p_face["radial_velocity"]) / 1e5 # km/s
    
    extent = (-width/2, width/2, -width/2, width/2)
    
    de_im = ax[1,0].imshow(d_edge, origin="lower", cmap='viridis',
                           norm=LogNorm(1e-32, 1e-24), extent=extent)
    te_im = ax[1,1].imshow(t_edge, origin="lower", cmap='magma',
                           norm=LogNorm(1e3,1e8), extent=extent)
    ve_im = ax[1,2].imshow(v_edge, origin="lower", cmap='coolwarm',
                           norm=SymLogNorm(1, linscale=0.5, base=10,
                                           vmin=-3e3, vmax=3e3),
                           extent=extent)
    df_im = ax[2,0].imshow(d_face, origin="lower", cmap='viridis',
                           norm=LogNorm(1e-32, 1e-24), extent=extent)
    tf_im = ax[2,1].imshow(t_face, origin="lower", cmap='magma',
                           norm=LogNorm(1e3,1e8), extent=extent)
    vf_im = ax[2,2].imshow(v_face, origin="lower", cmap='coolwarm',
                           norm=SymLogNorm(1, linscale=0.5, base=10,
                                           vmin=-3e3, vmax=3e3),
                           extent=extent)
    
    d_cb = plt.colorbar(de_im, cax=ax[0,0], orientation="horizontal")
    t_cb = plt.colorbar(te_im, cax=ax[0,1], orientation="horizontal")
    v_cb = plt.colorbar(ve_im, cax=ax[0,2], orientation="horizontal")
    
    d_cb.set_label(r'< Density >  [g cm$^{-3}$]', labelpad=-65)
    t_cb.set_label(r'< Temperature >  [K]', labelpad=-65)
    v_cb.set_label(r'< Radial Velocity >  [km/s]', labelpad=-65)
    
    for i in range(3):
        ax[0,i].tick_params(bottom=False, labelbottom=False,
                           top=True, labeltop=True)
    for i in range(1,3):
        ax[1,i].tick_params(labelleft=False, labelbottom=False)
        ax[2,i].tick_params(labelleft=False)
    ax[1,0].tick_params(labelbottom=False)
        
    #ax[1,0].set_xlabel("y  (kpc)")
    ax[1,0].set_ylabel("z  [kpc]")
    
    ax[2,0].set_xlabel("x  [kpc]")
    ax[2,0].set_ylabel("y  [kpc]")
    
    time = '{:.0f} Myr'.format(np.round(ds.current_time.to('Myr'),0))
    ax[2,0].text(0.03, 0.03, time, transform=ax[2,0].transAxes,
                 color='white', size='large', weight='bold')
  
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.05)
    fig.savefig(f'{ds.basename}_galaxy_ev.png')