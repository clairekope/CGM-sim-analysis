#!/usr/bin/env python
# coding: utf-8
############################################
# Use the function in plot_multi_profiles.py
# Converted from an ipython notebook
############################################

# In[1]:


from plot_multi_profiles import plot_profiles_over_time


# In[2]:


path = "../sim_images/lowres/mid_mass/profiles/"
title = 'mid_mass'


# In[3]:


start = 0
end = 1300
step = 50
dt = 50


# In[4]:


rad = 0
dens = 1
temp = 2
ent = 3
pres = 4
met = 5
tcool = 6
ratio = 7


# In[5]:


fig, ax = plot_profiles_over_time(rad, dens, start, end, step, 
                                  fpath=path, dtDD=dt, cmap='cividis', fbase='profiles')
ax.set_xlabel('r (kpc)')
ax.set_ylabel(r'$\rho$ (g/cc)')


# In[6]:


fig, ax = plot_profiles_over_time(rad, ent, start, end, step,
                                  fpath=path, dtDD=dt, fbase='profiles')
ax.set_xlabel('r (kpc)')
ax.set_ylabel('K (keV cm$^2$)')
# ax.axvline(1.3, c='k', ls='--') # 4*z_s
# ax.axvline(14, c='k', ls='--') # 4*r_s
ax.axvline(206, c='k', ls='-')
# ax.text(1.4, 5e2, r'$4z_s$')
# ax.text(15, 5e2, r'$4r_s$')
ax.text(216, 5e2, r'$r_\mathrm{vir}$')
# ax.set_xlim(1,500)
# ax.set_ylim(1e-3, 1e3)
# ax.annotate('Inital FB\nshockwave', (55, 80), (40, 0.8),
#             horizontalalignment='center',
#             arrowprops={'width':2,})

r = np.linspace(14,206)
K03 = 3.4 * r**0.71
K1 = 5.6 * r**0.71
voit1 = ax.plot(r,K03)
voit2 = ax.plot(r,K1)

ax.legend((voit1[0], voit2[0]),
          (r'$1.1\times10^{12}\ \mathrm{M_\odot}$,'+'\n'+r'$0.3\ \mathrm{Z_\odot}$ (Voit 19)',
           r'$1.1\times10^{12}\ \mathrm{M_\odot}$,'+'\n'+r'$1.0\ \mathrm{Z_\odot}$ (Voit 19)'))

ax.set_title(title)
fig.savefig(path+'entropy_evolution.png')


# In[7]:


fig, ax = plot_profiles_over_time(rad, tcool, start, end, step,
                                  fpath=path, dtDD=dt, fbase='profiles')
ax.set_xlabel('r (kpc)')
ax.set_ylabel(r'$t_\mathrm{cool}$ (Gyr)')
# ax.axvline(1.3, c='k', ls='--') # 4*z_s
# ax.axvline(14, c='k', ls='--') # 4*r_s
ax.axvline(206, c='k', ls='-')
# ax.text(1.4, 5e2, r'$4z_s$')
# ax.text(15, 5e2, r'$4r_s$')
ax.text(216, 5e2, r'$r_\mathrm{vir}$')
ax.axhline(1, c='k', ls=':')
# ax.set_xlim(1,500)
# ax.set_ylim(1e-3, 1e3)

ax.set_title(title)
fig.savefig(path+'tcool_evolution.png')


# In[8]:


fig, ax = plot_profiles_over_time(rad, ratio, start, end, step,
                                  fpath=path, dtDD=dt, fbase='profiles')
ax.set_xlabel('r (kpc)')
ax.set_ylabel(r'$t_\mathrm{cool}/t_\mathrm{ff}$')
# ax.axvline(1.3, c='k', ls='--') # 4*z_s
# ax.axvline(14, c='k', ls='--') # 4*r_s
ax.axvline(206, c='k', ls='-')
# ax.text(1.4, 7e2, r'$4z_s$')
# ax.text(15, 7e2, r'$4r_s$')
ax.text(216, 7e2, r'$r_\mathrm{vir}$')
ax.axhline(10, c='k', ls=':')
ax.axhline(20, c='k', ls=':')
# ax.set_xlim(1,500)
# ax.set_ylim(1e-1, 1e3)

ax.set_title(title)
fig.savefig(path+'tcool-tff_evolution.png')


# In[9]:


fig, ax = plot_profiles_over_time(rad, met, start, end, step,
                                  fpath=path, dtDD=dt, fbase='profiles')
ax.set_xlabel('r (kpc)')
ax.set_ylabel(r'$Z\ \mathrm{(Z_\odot)}}$')
# ax.axvline(1.3, c='k', ls='--') # 4*z_s
# ax.axvline(14, c='k', ls='--') # 4*r_s
ax.axvline(206, c='k', ls='-')
# ax.set_xlim(1,500)
# ax.set_ylim(1e-3, 1e3)


# In[ ]:




