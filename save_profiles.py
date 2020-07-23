import yt
yt.enable_parallelism()
import yt.units as u
import numpy as np

quants = ['density','temperature','entropy','pressure',
          'metallicity', 'cooling_time']

def NFW_mass_enclosed(radius):
    Mtot = 1e12 * u.Msun
    C = 10

    rho_crit = 1.8788e-29*0.49 * u.g/u.cm**3
    Rvir = ( 3.0/(4.0*np.pi)*Mtot / (200.*rho_crit) )**(1./3.)
    Rs = Rvir/C
    rho_0 = 200.0*rho_crit * C**3/3.0 / (np.log(1.0+C) - C/(1.0+C))

    M_r = 4.0*np.pi * rho_0 * Rs**3.0 * \
          (np.log((Rs+radius)/Rs) - radius/(Rs+radius))

    return M_r.to('g')

def MN_accel(radius):
    rs = 3.5 * u.kpc
    zs = 0.325 * u.kpc
    MStar = 5.8e10 * u.Msun
    G = u.G

    # sample some phi's and average them
    phi = np.linspace(0, np.pi/2) # symmetric
    phicol, radrow = np.meshgrid(phi,radius)
    r = radrow*np.sin(phicol) # cyl radius from sph
    z = radrow*np.cos(phicol)

    accel_r = G*MStar*r/np.power(np.power(r,2) \
              + np.power(
                  rs+np.sqrt(np.power(z,2) \
                             + zs**2),
                  2),3/2);

    accel_z = G*MStar*z/np.sqrt(np.power(z,2) \
              + zs**2)/np.power(np.power(r,2) \
              + np.power(rs+np.sqrt(np.power(z,2) \
                                    + zs**2),2),3/2) \
                *(  rs+np.sqrt(np.power(z,2) \
                  + zs**2)
                 );

    g = np.sqrt(np.power(accel_r,2) + np.power(accel_z,2)) # no phi accel

    return np.average(g, axis=1).to('cm/s**2')

if __name__=="__main__":
    datasets = yt.load("DD????/DD????")

    for ds in datasets.piter():

        sph = ds.sphere([0.5,0.5,0.5],(500,'kpc'))

        prof = yt.create_profile(sph, 'radius', quants, weight_field='cell_mass',
                                 units = {'radius':'kpc', 'cooling_time':'Gyr',
                                          'entropy':'keV*cm**2'})
        
        mass_enc = yt.create_profile(sph, 'radius', ['cell_mass'],
                                     weight_field=None, accumulation=True,
                                     units={'radius':'kpc',
                                            'cell_mass':'Msun'})

        assert np.allclose(prof.x.v, mass_enc.x.v)
        
        if ('all','particle_mass') in ds.field_list:
            mass_enc_part = yt.create_profile(sph, 'particle_radius', [('all','particle_mass')],
                                              weight_field=None, accumulation=True,
                                              units={'particle_radius':'kpc',
                                                     ('all','particle_mass'):'Msun'},
                                              override_bins={'particle_radius':prof.x_bins})

            assert np.allclose(prof.x.v, mass_enc_part.x.v)

        arr = np.empty((prof.x.size, len(quants)+2)) # add radius and tcool/tff
        arr[:,0] = prof.x
        for i, field in enumerate(quants):
            arr[:,i+1] = prof[field]

        # total NFW, gas, and particle mass. Average Miyamoto & Nagai acceleration
        Mtot_r = NFW_mass_enclosed(prof.x) + \
                 mass_enc['cell_mass']
        if ('all','particle_mass') in ds.field_list:
            Mtot_r += mass_enc_part[('all','particle_mass')]
        
        g_r = u.G * Mtot_r / np.power(prof.x,2)
        g_r = g_r.to('cm/s**2') + MN_accel(prof.x)
        
        arr[:,-1] = prof['cooling_time'] / np.sqrt(2*prof.x/g_r)

        
        header = 'radius'
        for field in quants:
            header += ' ' + field
        header += ' tcool/tff'

        np.savetxt('profiles_{:04d}.txt'.format(int(ds.basename[-4:])), arr,
                   header=header)
