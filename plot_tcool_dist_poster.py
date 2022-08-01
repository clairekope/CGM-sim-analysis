import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator

plt.rcParams.update({"font.size":18})

fid = np.genfromtxt("../extracted_data/fid_tcool_mass_dist_CGM.txt")
cflow = np.genfromtxt("../extracted_data/cflow_tcool_mass_dist_CGM.txt")
tctff5 = np.genfromtxt("../extracted_data/tctff5_tcool_mass_dist_CGM.txt")
tctff20 = np.genfromtxt("../extracted_data/tctff20_tcool_mass_dist_CGM.txt")
linrot = np.genfromtxt("../extracted_data/linrot_tcool_mass_dist_CGM.txt")
norot = np.genfromtxt("../extracted_data/norot_tcool_mass_dist_CGM.txt")

fid_disk = np.genfromtxt("../extracted_data/fid_tcool_mass_dist_CGM-disk.txt")
cflow_disk = np.genfromtxt("../extracted_data/cflow_tcool_mass_dist_CGM-disk.txt")
tctff5_disk = np.genfromtxt("../extracted_data/tctff5_tcool_mass_dist_CGM-disk.txt")
tctff20_disk = np.genfromtxt("../extracted_data/tctff20_tcool_mass_dist_CGM-disk.txt")
linrot_disk = np.genfromtxt("../extracted_data/linrot_tcool_mass_dist_CGM-disk.txt")
norot_disk = np.genfromtxt("../extracted_data/norot_tcool_mass_dist_CGM-disk.txt")

fig, ax = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(12,9))

ax[0,0].loglog(fid[:,0], np.cumsum(fid[:, 21]), c='C0', label='Fiducial')
ax[1,0].loglog(fid[:,0], np.cumsum(fid[:, 81]), c='C0', label='Fiducial')

ax[0,0].loglog(fid_disk[:,0], np.cumsum(fid_disk[:,1]), c='C0', ls='--')
ax[1,0].loglog(fid_disk[:,0], np.cumsum(fid_disk[:,5]), c='C0', ls='--')

ax[0,0].loglog(cflow[:,0], np.cumsum(cflow[:, 21]), c='C3', label='Cooling Flow')
ax[1,0].loglog(cflow[:,0], np.cumsum(cflow[:, 81]), c='C3', label='Cooling Flow')

ax[0,0].loglog(cflow_disk[:,0], np.cumsum(cflow_disk[:,1]), c='C3', ls='--')
ax[1,0].loglog(cflow_disk[:,0], np.cumsum(cflow_disk[:,5]), c='C3', ls='--')


ax[0,1].loglog(fid[:,0], np.cumsum(fid[:, 21]), c='C0')
ax[1,1].loglog(fid[:,0], np.cumsum(fid[:, 81]), c='C0')

ax[0,1].loglog(tctff5[:,0], np.cumsum(tctff5[:, 21]), c='C2', label=r'$t_{\rm c}/t_{\rm ff}=5$')
ax[1,1].loglog(tctff5[:,0], np.cumsum(tctff5[:, 81]), c='C2', label=r'$t_{\rm c}/t_{\rm ff}=5$')

ax[0,1].loglog(tctff5_disk[:,0], np.cumsum(tctff5_disk[:,1]), c='C2', ls='--')
ax[1,1].loglog(tctff5_disk[:,0], np.cumsum(tctff5_disk[:,5]), c='C2', ls='--')

ax[0,1].loglog(tctff20[:,0], np.cumsum(tctff20[:, 21]), c='C1', label=r'$t_{\rm c}/t_{\rm ff}=20$')
ax[1,1].loglog(tctff20[:,0], np.cumsum(tctff20[:, 81]), c='C1', label=r'$t_{\rm c}/t_{\rm ff}=20$')

ax[0,1].loglog(tctff20_disk[:,0], np.cumsum(tctff20_disk[:,1]), c='C1', ls='--')
ax[1,1].loglog(tctff20_disk[:,0], np.cumsum(tctff20_disk[:,5]), c='C1', ls='--')


ax[0,2].loglog(fid[:,0], np.cumsum(fid[:, 21]), c='C0')
ax[1,2].loglog(fid[:,0], np.cumsum(fid[:, 81]), c='C0')

ax[0,2].loglog(linrot[:,0], np.cumsum(linrot[:, 21]), c='C4', label='Linear Rot')
ax[1,2].loglog(linrot[:,0], np.cumsum(linrot[:, 81]), c='C4', label='Linear Rot')

ax[0,2].loglog(linrot_disk[:,0], np.cumsum(linrot_disk[:,1]), c='C4', ls='--')
ax[1,2].loglog(linrot_disk[:,0], np.cumsum(linrot_disk[:,5]), c='C4', ls='--')

ax[0,2].loglog(norot[:,0], np.cumsum(norot[:, 21]), c='C5', label='No Rot')
ax[1,2].loglog(norot[:,0], np.cumsum(norot[:, 81]), c='C5', label='No Rot')

ax[0,2].loglog(norot_disk[:,0], np.cumsum(norot_disk[:,1]), c='C5', ls='--')
ax[1,2].loglog(norot_disk[:,0], np.cumsum(norot_disk[:,5]), c='C5', ls='--')


#ax[0,0].set_title('0 Gyr', fontweight='bold')
#ax[1,0].set_title('4 Gyr', fontweight='bold')

ax[0,0].text(0.04, 0.88, "1 Gyr", transform=ax[0,0].transAxes,
             fontdict={"weight":"bold"})
ax[1,0].text(0.04, 0.88, "4 Gyr", transform=ax[1,0].transAxes,
             fontdict={"weight":"bold"})

cgm = Line2D([0], [0], ls='-', c='dimgray')
w_disk = Line2D([0], [0], ls='--', c='dimgray')

fid = Line2D([0], [0], ls='-', c='C0')
cflow = Line2D([0], [0], ls='-', c='C3')
tctff5 = Line2D([0], [0], ls='-', c='C2')
tctff20 = Line2D([0], [0], ls='-', c='C1')
linrot = Line2D([0], [0], ls='-', c='C4')
norot = Line2D([0], [0], ls='-', c='C5')

fig.legend([fid, cflow, tctff5, tctff20, linrot, norot, cgm, w_disk],
           ["Fiducial","Cooling Flow",r"$t_{\rm cool}/t_{\rm ff} = 5$",r"$t_{\rm cool}/t_{\rm ff} =20$",
            "Linear Rotation","No Rotation","CGM Only","CGM + Disk"],
           loc='upper center', ncol=3)

#ax[0,0].set_xlim(1e-5, 1e1)
ax[0,0].set_xlim(1e-1, 1e2)

for i in range(2):
    ax[i,0].set_ylabel(r"$M(t_{\rm c} < t_{\rm c,0})$", fontsize='large')
    ax[i,0].yaxis.set_minor_locator(FixedLocator([1e8,1e10]))

for i in range(3):
    ax[1,i].set_xlabel(r"$t_{\rm c,0}$  [Gyr]", fontsize='large')
        
ax[0,0].set_ylim(1e7, 1e11)

for i in range(2):
    for j in range(3):
        ax[i,j].grid(axis='both')
        ax[i,j].grid(axis='y', which='minor')

fig.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.82)
fig.savefig("../fig_tcool-mass-dist_cumm-big_poster.pdf", transparent=True)
