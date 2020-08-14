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

def cell_and_particle_mass_enclosed(data_obj, rad_bins, rad_unit='kpc'):
    """
    Returns cummulative cell (and particle) mass from 
    yt data container data_obj within radial bins rad_bins

    Example usage:

    sph = ds.sphere([0.5,0.5,0.5], (500,'kpc'))
    prof = yt.create_profile(sph, 'radius', 'density', units={'radius':'kpc'})
    M_enc = gas_and_particle_mass_enclosed(sph, prof.x_bins.v)
    """

    mass_enc = yt.create_profile(data_obj, 'radius', ['cell_mass'],
                                 weight_field=None, accumulation=True,
                                 units={'radius':'kpc',
                                        'cell_mass':'Msun'},
                                 override_bins={'radius':rad_bins})

    assert np.allclose(rad_bins, mass_enc.x_bins.v)
        
    if ('all','particle_mass') in data_obj.ds.field_list:
        mass_enc_part = yt.create_profile(data_obj,
                                          'particle_radius', [('all','particle_mass')],
                                          weight_field=None, accumulation=True,
                                          units={'particle_radius':'kpc',
                                                 ('all','particle_mass'):'Msun'},
                                          override_bins={'particle_radius':rad_bins})

        assert np.allclose(rad_bins, mass_enc_part.x_bins.v)

    M_enc = mass_enc['cell_mass']
    if ('all','particle_mass') in data_obj.ds.field_list:
        M_enc += mass_enc_part[('all','particle_mass')]

    return M_enc


if __name__=="__main__":
    datasets = yt.load("DD????/DD????")

    if os.path.isfile('last_profile'):
        with open('last_profile','r') as f:
            last_profile = f.read()

    else:
        last_profile = -1

    for ds in datasets.piter():
        if int(ds.basename[-4:]) > last_profile:
            center = ds.quan(0.5, 'code_length')
            le = center-ds.quan(0.6, 'Mpc')
            re = center+ds.quan(0.6, 'Mpc')
            offset = ds.quan(5, 'kpc')

            sph = ds.sphere('c',(500,'kpc'))        
            dsk_slb = ds.region([center, center, center],
                                [le, le, center-offset],
                                [re, re, center+offset])

            cgm = sph-dsk_slb
            cgm.set_field_parameter('center', ds.arr([0.5,0.5,0.5],'code_length'))

            prof = yt.create_profile(cgm, 'radius', quants, weight_field='cell_mass',
                                     units = {'radius':'kpc', 'cooling_time':'Gyr',
                                              'entropy':'keV*cm**2'})

            arr = np.empty((prof.x.size, len(quants)+2)) # add radius and tcool/tff
            arr[:,0] = prof.x
            for i, field in enumerate(quants):
                arr[:,i+1] = prof[field]

            # total NFW, gas, and particle mass. Average Miyamoto & Nagai acceleration
            Mtot_r = NFW_mass_enclosed(prof.x) + \
                     cell_and_particle_mass_enclosed(cgm, prof.x_bins.v)
        
            g_r = u.G * Mtot_r / np.power(prof.x,2)
            g_r = g_r.to('cm/s**2') + MN_accel(prof.x)
        
            arr[:,-1] = prof['cooling_time'] / np.sqrt(2*prof.x/g_r)
        
            header = 'radius'
            for field in quants:
                header += ' ' + field
                header += ' tcool/tff'

            np.savetxt('profiles_{:04d}.txt'.format(int(ds.basename[-4:])), arr,
                       header=header)

    last_profile = int(ds[-1].basename[-4:])
    with open('last_profile','w') as f:
        f.write(last_profile)
