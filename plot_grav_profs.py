# coding: utf-8
import unyt as u
import numpy as np
import matplotlib.pyplot as plt
from calc_toomre import DiskModel

r = np.linspace(1,201) * u.kpc

dsk = DiskModel()

g_nfw_mn = dsk.g(r)
g_modnfw = dsk.g_modNFW(r)

fig, ax = plt.subplots(ncols=2)

ax[0].plot(r, np.sqrt(r*g_nfw_mn).to('km/s'), label='NFW+MN')
ax[0].plot(r, np.sqrt(r*g_modnfw).to('km/s'), label='modNFW')

ax[1].plot(r, g_nfw_mn.to('km/s**2'))
ax[1].plot(r, g_modnfw.to('km/s**2'))

ax[0].set_xlabel("r [kpc]")
ax[1].set_xlabel("r [kpc]")

ax[0].set_ylabel(r"$v_\mathrm{circ}$ [km/s]")
ax[1].set_ylabel(r"g [$\mathrm{km/s^2}$]")

ax[0].legend()
plt.show()
