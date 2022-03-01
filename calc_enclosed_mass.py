#####################################################
# Calculate mass enclosed as a function of radius for
# 1. NFW dark matter profile
# 2. Miyamoto & Nagai stellar profile (angle averaged)
# 3. Simulation cells & particles
#####################################################

import yt
import yt.units as u
import numpy as np

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

    # sample some theta's and average them
    theta = np.linspace(0, np.pi/2) # symmetric
    thetacol, radrow = np.meshgrid(theta,radius)
    r = radrow*np.sin(thetacol) # cyl radius from sph
    z = radrow*np.cos(thetacol)

    accel_r = G*MStar*r/np.power(np.power(r,2) \
              + np.power(
                  rs+np.sqrt(np.power(z,2) \
                             + zs**2),
                  2),3/2)

    accel_z = G*MStar*z/np.sqrt(np.power(z,2) \
              + zs**2)/np.power(np.power(r,2) \
              + np.power(rs+np.sqrt(np.power(z,2) \
                                    + zs**2),2),3/2) \
                *(  rs+np.sqrt(np.power(z,2) \
                  + zs**2)
                 )

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
                                 units={'radius':rad_unit,
                                        'cell_mass':'Msun'},
                                 override_bins={'radius':rad_bins})

    assert np.allclose(rad_bins, mass_enc.x_bins.v)
    
    M_enc = mass_enc['cell_mass']
        
    if ('all','particle_mass') in data_obj.ds.field_list:
        mass_enc_part = yt.create_profile(data_obj,
                                          'particle_radius', [('all','particle_mass')],
                                          weight_field=None, accumulation=True,
                                          units={'particle_radius':rad_unit,
                                                 ('all','particle_mass'):'Msun'},
                                          override_bins={'particle_radius':rad_bins})

        assert np.allclose(rad_bins, mass_enc_part.x_bins.v)

        M_enc += mass_enc_part[('all','particle_mass')]

    return M_enc.to('g')

