import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt

out = "DD0000"
ds = out+"/"+out

ds_md = yt.load("GS_lowres_fid/"+ds)
ds_lo = yt.load("GS_lowres_low_mass/"+ds)
ds_hi = yt.load("GS_lowres_high_mass/"+ds)

sph_md = ds_md.sphere('c', (500, 'kpc'))
sph_lo = ds_lo.sphere('c', (500, 'kpc'))
sph_hi = ds_hi.sphere('c', (500, 'kpc'))

prof_md = yt.create_profile(sph_md, 'radius', ['density','temperature','entropy','pressure'],
                            units={'radius':'kpc'}, logs={'radius':False})
prof_lo = yt.create_profile(sph_lo, 'radius', ['density','temperature','entropy','pressure'],
                            units={'radius':'kpc'}, logs={'radius':False})
prof_hi = yt.create_profile(sph_hi, 'radius', ['density','temperature','entropy','pressure'],
                            units={'radius':'kpc'}, logs={'radius':False})

fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(7,9))

ax[0].set_ylabel(r"$\mathrm{\rho\ [g\ cm^{-3}]}$")
ax[1].set_ylabel(r"T [K]")
ax[2].set_ylabel(r"$\mathrm{K\ [cm^2 \cdot keV]}$")

#ax[0].set_xlabel("r [kpc]")
#ax[1].set_xlabel("r [kpc]")
ax[2].set_xlabel("r [kpc]")

ax[0].loglog(prof_md.x, prof_md['density'], label=r"$\mathrm{M_{fid}}$")
ax[0].loglog(prof_lo.x, prof_lo['density'], ":", label=r"$\mathrm{M_{fid}}/3$")
ax[0].loglog(prof_hi.x, prof_hi['density'], "--", label=r"$\mathrm{3M_{fid}}$")

ax[1].loglog(prof_md.x, prof_md['temperature'], label=r"$\mathrm{M_{fid}}")
ax[1].loglog(prof_lo.x, prof_lo['temperature'], ":", label=r"$\mathrm{M_{fid}}/3$")
ax[1].loglog(prof_hi.x, prof_hi['temperature'], "--", label=r"$\mathrm{3M_{fid}}$")

ax[2].loglog(prof_md.x, prof_md['entropy'], label=r"$\mathrm{M_{fid}}$")
ax[2].loglog(prof_lo.x, prof_lo['entropy'], ":", label=r"$\mathrm{M_{fid}}/3$")
ax[2].loglog(prof_hi.x, prof_hi['entropy'], "--", label=r"$\mathrm{3M_{fid}}$")

ax[3].loglog(prof_md.x, prof_md['pressure'], label=r"$\mathrm{M_{fid}}$")
ax[3].loglog(prof_lo.x, prof_lo['pressure'], ":", label=r"$\mathrm{M_{fid}}/3$")
ax[3].loglog(prof_hi.x, prof_hi['pressure'], "--", label=r"$\mathrm{3M_{fid}}$")

ax[0].legend()

fig.tight_layout()
fig.savefig("mass_variant_profiles_{:s}".format(ds_md.basename))
