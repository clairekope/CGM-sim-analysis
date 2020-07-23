import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# velocity in km/s, length in kpc

def vmag(x, y, v0=180, r0=10, beta=-0.5):
    r = np.sqrt(x**2+y**2)
    theta = np.arccos(y/r) # z/r
    v = np.sin(theta)**2 \
        * np.where(r <= r0, v0, v0*(r/r0)**beta)
    return v

X = np.linspace(2, 200, 500)
x, y = np.meshgrid(X,X)

p = plt.pcolormesh(X, X, vmag(x,y), norm=LogNorm(), vmax=200, vmin=0.1)
c = plt.contour(X, X, vmag(x,y), 10, cmap='bone')
cb = plt.colorbar(p)

plt.clabel(c, inline=1, colors='k', fontsize=10, fmt="%1.0f")
cb.set_label('km/s')

plt.xlabel('$\mathrm{r_{cyl}}$ (kpc)')
plt.ylabel('z (kpc)')
plt.show()
