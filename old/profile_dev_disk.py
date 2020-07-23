import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

TruncRadius = 31.2 * 3e21 # kpc to cm
R0 = 3.5 * 3e21 # kpc to cm
Z0 = 0.325 * 3e21 # kpc to cm
Rho0 = 4.9e9 * 1.988e33 / (8*np.pi*Z0*R0**2) # Msun/cc to g/cc

SmoothRadius = TruncRadius*.02/.026
SmoothLength = TruncRadius - SmoothRadius

def TaskerDensity(r, z, rho0=Rho0, r0=R0, z0=Z0):
    return rho0 * np.exp(-r/r0) / np.cosh(0.5*z/z0)**2 

def SechDensity(r, z, rho0=Rho0, r0=R0, z0=Z0):
    dens = np.where(r < TruncRadius, 0.25 * rho0 / np.cosh(r/r0) / np.cosh(z/z0), 0.0)
    return dens

def TonnesenDensity(r, z, rho0=Rho0, r0=R0, z0=Z0):
    dens = SechDensity(r, z, rho0, r0, z0)
    dens = np.where(np.logical_and(SmoothRadius < r, r < TruncRadius),
                    dens * 0.5*(1+np.cos(np.pi * (r-SmoothRadius)/SmoothLength)),
                    dens)

    return dens

def DoubleExpDensity(r, z, rho0=Rho0, r0=R0, z0=Z0):
    return rho0 * np.exp(-r/r0) * np.exp(-np.abs(z)/z0)

X = np.linspace(0, 35) # kpc
Y = np.linspace(-10, 10) # kpc
x, y = np.meshgrid(X*3e21, Y*3e21) # cm
r = np.sqrt(x**2 + y**2) # cm

fig, ax = plt.subplots(nrows=3, sharex=True, sharey=True)
lvls = np.linspace(-9,1,10+1)

c = ax[0].contourf(X, Y, np.log10(DoubleExpDensity(r,y)/1.67e-24), levels=lvls)
ax[0].set_title("Double Exponential")

ax[1].contourf(X, Y, np.log10(TaskerDensity( r,y)/1.67e-24), levels=lvls)
ax[1].set_title("Tasker & Bryan '06")

ax[2].contourf(X, Y, np.log10(TonnesenDensity(r,y)/1.67e-24), levels=lvls)
ax[2].set_title("Tonnesen & Bryan '09")

ax[0].set_ylabel('z (kpc)')
ax[1].set_ylabel('z (kpc)')
ax[2].set_ylabel('z (kpc)')
ax[2].set_xlabel('r (kpc)')
fig.tight_layout()

cb = plt.colorbar(c, ax=ax)
cb.set_label('log(particles/cc)')

plt.show()
