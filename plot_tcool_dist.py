import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

fid = np.genfromtxt("../original_sims/extracted_data/fid_tcool_mass_dist_CGM.txt")
cflow = np.genfromtxt("../original_sims/extracted_data/cflow_tcool_mass_dist_CGM.txt")
tctff5 = np.genfromtxt("../original_sims/extracted_data/tctff5_tcool_mass_dist_CGM.txt")
tctff20 = np.genfromtxt("../original_sims/extracted_data/tctff20_tcool_mass_dist_CGM.txt")
linrot = np.genfromtxt("../original_sims/extracted_data/linrot_tcool_mass_dist_CGM.txt")
norot = np.genfromtxt("../original_sims/extracted_data/norot_tcool_mass_dist_CGM.txt")

fid_disk = np.genfromtxt("../original_sims/extracted_data/fid_tcool_mass_dist_CGM-disk.txt")
cflow_disk = np.genfromtxt("../original_sims/extracted_data/cflow_tcool_mass_dist_CGM-disk.txt")
tctff5_disk = np.genfromtxt("../original_sims/extracted_data/tctff5_tcool_mass_dist_CGM-disk.txt")
tctff20_disk = np.genfromtxt("../original_sims/extracted_data/tctff20_tcool_mass_dist_CGM-disk.txt")
linrot_disk = np.genfromtxt("../original_sims/extracted_data/linrot_tcool_mass_dist_CGM-disk.txt")
norot_disk = np.genfromtxt("../original_sims/extracted_data/norot_tcool_mass_dist_CGM-disk.txt")

fig, ax = plt.subplots(nrows=3, ncols=5, sharex=True, sharey=True, figsize=(11,7))

ax[0,0].loglog(fid[:,0], np.cumsum(fid[:,  1]), c='C0', label='Fiducial')
ax[0,1].loglog(fid[:,0], np.cumsum(fid[:, 21]), c='C0')
ax[0,2].loglog(fid[:,0], np.cumsum(fid[:, 41]), c='C0')
ax[0,3].loglog(fid[:,0], np.cumsum(fid[:, 61]), c='C0')
ax[0,4].loglog(fid[:,0], np.cumsum(fid[:, 81]), c='C0', label='Fiducial')

ax[0,0].loglog(fid_disk[:,0], np.cumsum(fid_disk[:,1]), c='C0', ls='--')
ax[0,1].loglog(fid_disk[:,0], np.cumsum(fid_disk[:,2]), c='C0', ls='--')
ax[0,2].loglog(fid_disk[:,0], np.cumsum(fid_disk[:,3]), c='C0', ls='--')
ax[0,3].loglog(fid_disk[:,0], np.cumsum(fid_disk[:,4]), c='C0', ls='--')
ax[0,4].loglog(fid_disk[:,0], np.cumsum(fid_disk[:,5]), c='C0', ls='--')

ax[0,0].loglog(cflow[:,0], np.cumsum(cflow[:,  1]), c='C3', label='Cooling Flow')
ax[0,1].loglog(cflow[:,0], np.cumsum(cflow[:, 21]), c='C3')
ax[0,2].loglog(cflow[:,0], np.cumsum(cflow[:, 41]), c='C3')
ax[0,3].loglog(cflow[:,0], np.cumsum(cflow[:, 61]), c='C3')
ax[0,4].loglog(cflow[:,0], np.cumsum(cflow[:, 81]), c='C3', label='Cooling Flow')

ax[0,0].loglog(cflow_disk[:,0], np.cumsum(cflow_disk[:,1]), c='C3', ls='--')
ax[0,1].loglog(cflow_disk[:,0], np.cumsum(cflow_disk[:,2]), c='C3', ls='--')
ax[0,2].loglog(cflow_disk[:,0], np.cumsum(cflow_disk[:,3]), c='C3', ls='--')
ax[0,3].loglog(cflow_disk[:,0], np.cumsum(cflow_disk[:,4]), c='C3', ls='--')
ax[0,4].loglog(cflow_disk[:,0], np.cumsum(cflow_disk[:,5]), c='C3', ls='--')

ax[0,0].legend(loc='lower right', framealpha=1)

ax[1,0].loglog(fid[:,0], np.cumsum(fid[:,  1]), c='C0')
ax[1,1].loglog(fid[:,0], np.cumsum(fid[:, 21]), c='C0')
ax[1,2].loglog(fid[:,0], np.cumsum(fid[:, 41]), c='C0')
ax[1,3].loglog(fid[:,0], np.cumsum(fid[:, 61]), c='C0')
ax[1,4].loglog(fid[:,0], np.cumsum(fid[:, 81]), c='C0')

ax[1,0].loglog(tctff5[:,0], np.cumsum(tctff5[:,  1]), c='C2', label=r'$t_{\rm c}/t_{\rm ff}=5$')
ax[1,1].loglog(tctff5[:,0], np.cumsum(tctff5[:, 21]), c='C2')
ax[1,2].loglog(tctff5[:,0], np.cumsum(tctff5[:, 41]), c='C2')
ax[1,3].loglog(tctff5[:,0], np.cumsum(tctff5[:, 61]), c='C2')
ax[1,4].loglog(tctff5[:,0], np.cumsum(tctff5[:, 81]), c='C2', label=r'$t_{\rm c}/t_{\rm ff}=5$')

ax[1,0].loglog(tctff5_disk[:,0], np.cumsum(tctff5_disk[:,1]), c='C2', ls='--')
ax[1,1].loglog(tctff5_disk[:,0], np.cumsum(tctff5_disk[:,2]), c='C2', ls='--')
ax[1,2].loglog(tctff5_disk[:,0], np.cumsum(tctff5_disk[:,3]), c='C2', ls='--')
ax[1,3].loglog(tctff5_disk[:,0], np.cumsum(tctff5_disk[:,4]), c='C2', ls='--')
ax[1,4].loglog(tctff5_disk[:,0], np.cumsum(tctff5_disk[:,5]), c='C2', ls='--')

ax[1,0].loglog(tctff20[:,0], np.cumsum(tctff20[:,  1]), c='C1', label=r'$t_{\rm c}/t_{\rm ff}=20$')
ax[1,1].loglog(tctff20[:,0], np.cumsum(tctff20[:, 21]), c='C1')
ax[1,2].loglog(tctff20[:,0], np.cumsum(tctff20[:, 41]), c='C1')
ax[1,3].loglog(tctff20[:,0], np.cumsum(tctff20[:, 61]), c='C1')
ax[1,4].loglog(tctff20[:,0], np.cumsum(tctff20[:, 81]), c='C1', label=r'$t_{\rm c}/t_{\rm ff}=20$')

ax[1,0].loglog(tctff20_disk[:,0], np.cumsum(tctff20_disk[:,1]), c='C1', ls='--')
ax[1,1].loglog(tctff20_disk[:,0], np.cumsum(tctff20_disk[:,2]), c='C1', ls='--')
ax[1,2].loglog(tctff20_disk[:,0], np.cumsum(tctff20_disk[:,3]), c='C1', ls='--')
ax[1,3].loglog(tctff20_disk[:,0], np.cumsum(tctff20_disk[:,4]), c='C1', ls='--')
ax[1,4].loglog(tctff20_disk[:,0], np.cumsum(tctff20_disk[:,5]), c='C1', ls='--')

ax[1,0].legend(loc='lower right', framealpha=1)

ax[2,0].loglog(fid[:,0], np.cumsum(fid[:,  1]), c='C0')
ax[2,1].loglog(fid[:,0], np.cumsum(fid[:, 21]), c='C0')
ax[2,2].loglog(fid[:,0], np.cumsum(fid[:, 41]), c='C0')
ax[2,3].loglog(fid[:,0], np.cumsum(fid[:, 61]), c='C0')
ax[2,4].loglog(fid[:,0], np.cumsum(fid[:, 81]), c='C0')

ax[2,0].loglog(linrot[:,0], np.cumsum(linrot[:,  1]), c='C4', label='Linear Rot')
ax[2,1].loglog(linrot[:,0], np.cumsum(linrot[:, 21]), c='C4')
ax[2,2].loglog(linrot[:,0], np.cumsum(linrot[:, 41]), c='C4')
ax[2,3].loglog(linrot[:,0], np.cumsum(linrot[:, 61]), c='C4')
ax[2,4].loglog(linrot[:,0], np.cumsum(linrot[:, 81]), c='C4', label='Linear Rot')

ax[2,0].loglog(linrot_disk[:,0], np.cumsum(linrot_disk[:,1]), c='C4', ls='--')
ax[2,1].loglog(linrot_disk[:,0], np.cumsum(linrot_disk[:,2]), c='C4', ls='--')
ax[2,2].loglog(linrot_disk[:,0], np.cumsum(linrot_disk[:,3]), c='C4', ls='--')
ax[2,3].loglog(linrot_disk[:,0], np.cumsum(linrot_disk[:,4]), c='C4', ls='--')
ax[2,4].loglog(linrot_disk[:,0], np.cumsum(linrot_disk[:,5]), c='C4', ls='--')

ax[2,0].loglog(norot[:,0], np.cumsum(norot[:,  1]), c='C5', label='No Rot')
ax[2,1].loglog(norot[:,0], np.cumsum(norot[:, 21]), c='C5')
ax[2,2].loglog(norot[:,0], np.cumsum(norot[:, 41]), c='C5')
ax[2,3].loglog(norot[:,0], np.cumsum(norot[:, 61]), c='C5')
ax[2,4].loglog(norot[:,0], np.cumsum(norot[:, 81]), c='C5', label='No Rot')

ax[2,0].loglog(norot_disk[:,0], np.cumsum(norot_disk[:,1]), c='C5', ls='--')
ax[2,1].loglog(norot_disk[:,0], np.cumsum(norot_disk[:,2]), c='C5', ls='--')
ax[2,2].loglog(norot_disk[:,0], np.cumsum(norot_disk[:,3]), c='C5', ls='--')
ax[2,3].loglog(norot_disk[:,0], np.cumsum(norot_disk[:,4]), c='C5', ls='--')
ax[2,4].loglog(norot_disk[:,0], np.cumsum(norot_disk[:,5]), c='C5', ls='--')

ax[2,0].legend(loc='lower right', framealpha=1)

ax[0,0].set_title('0 Gyr', fontweight='bold')
ax[0,1].set_title('1 Gyr', fontweight='bold')
ax[0,2].set_title('2 Gyr', fontweight='bold')
ax[0,3].set_title('3 Gyr', fontweight='bold')
ax[0,4].set_title('4 Gyr', fontweight='bold')

cgm = Line2D([0], [0], ls='-', c='dimgray')
w_disk = Line2D([0], [0], ls='--', c='dimgray')

ax[0,1].legend([cgm, w_disk], ["CGM Only","CGM + Disk"])

#ax[0,0].set_xlim(1e-5, 1e1)
ax[0,0].set_xlim(1e-1, 1e2)

for i in range(3):
    ax[i,0].set_ylabel(r"$M(t_{\rm c} < t_{\rm c,0})\ \mathrm{[M_\odot]}$", fontsize='large')
#     for j in range(5):
#         ax[i,j].grid()
#         ax[i,j].axvline(1, c='gray', ls='--')
#         ax[i,j].axvline(10, c='gray', ls='--')
for i in range(5):
    ax[2,i].set_xlabel(r"$t_{\rm c,0}$  [Gyr]", fontsize='large')
        
ax[0,0].set_ylim(1e7, 1e11)
#ax[0,0].yaxis.set_major_locator(FixedLocator([1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10]))
#ax[0,0].set_xticks([1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])
#ax[0,0].set_xticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5])
#ax[0,0].set_xticklabels(['$10^{-2}$','','$10^{0}$','','$10^{2}$','','$10^4$',''])

for i in range(3):
    for j in range(5):
        ax[i,j].grid(axis='both')

fig.tight_layout()
fig.savefig("../original_sims/figures/fig_tcool-mass-dist_cumm-big.pdf")



fig, ax = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(5,5))

lb = 26
ub = -20

ax[0,0].loglog(fid[:,0], 
               np.cumsum(fid[:,  1])/np.cumsum(fid[:,  1])[-1],
               c='C0', label='Fid')
ax[0,1].loglog(fid[:,0], 
               np.cumsum(fid[:, 21])/np.cumsum(fid[:, 21])[-1],
               c='C0', label='Fid')
ax[1,0].loglog(fid[:,0], 
               np.cumsum(fid[:, 41])/np.cumsum(fid[:, 41])[-1],
               c='C0')
ax[1,1].loglog(fid[:,0], 
               np.cumsum(fid[:, 61])/np.cumsum(fid[:, 61])[-1],
               c='C0')
ax[2,0].loglog(fid[:,0], 
               np.cumsum(fid[:, 81])/np.cumsum(fid[:, 81])[-1],
               c='C0')

ax[0,0].loglog(tctff5[:,0],
               np.cumsum(tctff5[:,  1]) / np.cumsum(tctff5[:,  1])[-1],
               c='C2', label='LowRatio')
ax[0,1].loglog(tctff5[:,0],
               np.cumsum(tctff5[:, 21]) / np.cumsum(tctff5[:, 21])[-1],
               c='C2', label='LowRatio')
ax[1,0].loglog(tctff5[:,0],
               np.cumsum(tctff5[:, 41]) / np.cumsum(tctff5[:, 41])[-1],
               c='C2')
ax[1,1].loglog(tctff5[:,0],
               np.cumsum(tctff5[:, 61]) / np.cumsum(tctff5[:, 61])[-1],
               c='C2')
ax[2,0].loglog(tctff5[:,0],
               np.cumsum(tctff5[:, 81]) / np.cumsum(tctff5[:, 81])[-1],
               c='C2')

ax[0,0].loglog(tctff20[:,0],
               np.cumsum(tctff20[:,  1]) / np.cumsum(tctff20[:,  1])[-1],
               c='C1', label='HighRatio')
ax[0,1].loglog(tctff20[:,0],
               np.cumsum(tctff20[:, 21]) / np.cumsum(tctff20[:, 21])[-1],
               c='C1', label='HighRatio')
ax[1,0].loglog(tctff20[:,0],
               np.cumsum(tctff20[:, 41]) / np.cumsum(tctff20[:, 41])[-1],
               c='C1')
ax[1,1].loglog(tctff20[:,0],
               np.cumsum(tctff20[:, 61]) / np.cumsum(tctff20[:, 61])[-1],
               c='C1')
ax[2,0].loglog(tctff20[:,0],
               np.cumsum(tctff20[:, 81]) / np.cumsum(tctff20[:, 81])[-1],
               c='C1')

ax[0,1].legend(framealpha=1)

ax[0,0].text(0.05, 0.9, '0 Gyr', transform=ax[0,0].transAxes, fontweight='bold',
             ha='left', va='top')
ax[0,1].text(0.05, 0.9, '1 Gyr', transform=ax[0,1].transAxes, fontweight='bold',
             ha='left', va='top')
ax[1,0].text(0.05, 0.9, '2 Gyr', transform=ax[1,0].transAxes, fontweight='bold',
             ha='left', va='top')
ax[1,1].text(0.05, 0.9, '3 Gyr', transform=ax[1,1].transAxes, fontweight='bold',
             ha='left', va='top')
ax[2,0].text(0.05, 0.9, '4 Gyr', transform=ax[2,0].transAxes, fontweight='bold',
             ha='left', va='top')

ax[0,0].set_xlim(1e-1, 1e2)
ax[0,0].set_ylim(1e-3,1.5e0)

for i in range(2):
    ax[i,0].set_ylabel(r"Mass PDF", fontsize='large')
    for j in range(2):
        ax[i,j].grid()
        
ax[2,0].grid()
ax[2,1].set_axis_off()
ax[1,1].tick_params(labelbottom=True)

ax[2,0].set_ylabel(r"Mass PDF", fontsize='large')
ax[2,0].set_xlabel(r"$t_{\rm c,0}$  [Gyr]", fontsize='large')
ax[1,1].set_xlabel(r"$t_{\rm c,0}$  [Gyr]", fontsize='large')

fig.tight_layout()
fig.subplots_adjust(hspace=0.1)
fig.savefig("../original_sims/figures/fig_tcool-mass-dist_tctff-var.pdf")



fig, ax = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(5,5))

lb = 26
ub = -20

ax[0,0].loglog(fid[:,0], 
               np.cumsum(fid[:,  1])/np.cumsum(fid[:,  1])[-1],
               c='C0', label='Fid')
ax[0,1].loglog(fid[:,0], 
               np.cumsum(fid[:, 21])/np.cumsum(fid[:, 21])[-1],
               c='C0', label='Fid')
ax[1,0].loglog(fid[:,0], 
               np.cumsum(fid[:, 41])/np.cumsum(fid[:, 41])[-1],
               c='C0')
ax[1,1].loglog(fid[:,0], 
               np.cumsum(fid[:, 61])/np.cumsum(fid[:, 61])[-1],
               c='C0')
ax[2,0].loglog(fid[:,0], 
               np.cumsum(fid[:, 81])/np.cumsum(fid[:, 81])[-1],
               c='C0')

ax[0,0].loglog(linrot[:,0],
               np.cumsum(linrot[:,  1]) / np.cumsum(linrot[:,  1])[-1],
               c='C4', label='LinRot')
ax[0,1].loglog(linrot[:,0],
               np.cumsum(linrot[:, 21]) / np.cumsum(linrot[:, 21])[-1],
               c='C4', label='LinRot')
ax[1,0].loglog(linrot[:,0],
               np.cumsum(linrot[:, 41]) / np.cumsum(linrot[:, 41])[-1],
               c='C4')
ax[1,1].loglog(linrot[:,0],
               np.cumsum(linrot[:, 61]) / np.cumsum(linrot[:, 61])[-1],
               c='C4')
ax[2,0].loglog(linrot[:,0],
               np.cumsum(linrot[:, 81]) / np.cumsum(linrot[:, 81])[-1],
               c='C4')

ax[0,0].loglog(norot[:,0],
               np.cumsum(norot[:,  1]) / np.cumsum(norot[:,  1])[-1],
               c='C5', label='NoRot')
ax[0,1].loglog(norot[:,0],
               np.cumsum(norot[:, 21]) / np.cumsum(norot[:, 21])[-1],
               c='C5', label='NoRot')
ax[1,0].loglog(norot[:,0],
               np.cumsum(norot[:, 41]) / np.cumsum(norot[:, 41])[-1],
               c='C5')
ax[1,1].loglog(norot[:,0],
               np.cumsum(norot[:, 61]) / np.cumsum(norot[:, 61])[-1],
               c='C5')
ax[2,0].loglog(norot[:,0],
               np.cumsum(norot[:, 81]) / np.cumsum(norot[:, 81])[-1],
               c='C5')

ax[0,1].legend(framealpha=1)

ax[0,0].text(0.05, 0.9, '0 Gyr', transform=ax[0,0].transAxes, fontweight='bold',
             ha='left', va='top')
ax[0,1].text(0.05, 0.9, '1 Gyr', transform=ax[0,1].transAxes, fontweight='bold',
             ha='left', va='top')
ax[1,0].text(0.05, 0.9, '2 Gyr', transform=ax[1,0].transAxes, fontweight='bold',
             ha='left', va='top')
ax[1,1].text(0.05, 0.9, '3 Gyr', transform=ax[1,1].transAxes, fontweight='bold',
             ha='left', va='top')
ax[2,0].text(0.05, 0.9, '4 Gyr', transform=ax[2,0].transAxes, fontweight='bold',
             ha='left', va='top')

ax[0,0].set_xlim(1e-1, 1e2)
ax[0,0].set_ylim(1e-3,1.5e0)

for i in range(2):
    ax[i,0].set_ylabel(r"Mass PDF", fontsize='large')
    for j in range(2):
        ax[i,j].grid()
        
ax[2,0].grid()
ax[2,1].set_axis_off()
ax[1,1].tick_params(labelbottom=True)

ax[2,0].set_ylabel(r"Mass PDF", fontsize='large')
ax[2,0].set_xlabel(r"$t_{\rm c,0}$  [Gyr]", fontsize='large')
ax[1,1].set_xlabel(r"$t_{\rm c,0}$  [Gyr]", fontsize='large')

fig.tight_layout()
fig.subplots_adjust(hspace=0.1)
fig.savefig("../original_sims/figures/fig_tcool-mass-dist_rot-var.pdf")
