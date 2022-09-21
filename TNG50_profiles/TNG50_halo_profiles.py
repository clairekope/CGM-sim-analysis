# I'm reusing code from my breakBRD project but I don't want to store it in the
# same repository, so I have copied over functions from utilities.py instead of
# importing them. This script itself is modified from green_valley_properties.py
# I will also be adding duplicates of readsubfHDF5, readhaloHDF5, & snapHDF5
# to the repository for this purpose

import gc
import warnings
import requests
import readsubfHDF5
import readhaloHDF5
import snapHDF5
import numpy as np
import astropy.units as u
from astropy.constants import m_p, k_B, G
from mpi4py import MPI

warnings.simplefilter(action='ignore', category=FutureWarning)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


#
# Define utility functions
#

def scatter_work(array, mpi_rank, mpi_size, root=0, dtype=np.int32):
    """ will only work if MPI has been initialized by calling script.
        array should only exist on root & be None elsewhere"""
    if mpi_rank == root:
        scatter_total = array.size
        mod = scatter_total % mpi_size
        if mod != 0:
            print("Padding array for scattering...")
            pad = -1 * np.ones(mpi_size - mod, dtype=dtype)
            array = np.concatenate((array, pad))
            scatter_total += mpi_size - mod
            assert scatter_total % mpi_size == 0
            assert scatter_total == array.size
    else:
        scatter_total = None

    scatter_total = comm.bcast(scatter_total, root=root)
    subset = np.empty(scatter_total//mpi_size, dtype=dtype)
    comm.Scatter(array, subset, root=root)

    return subset

def get(path, params=None, fpath=""):
    
    attempt = 1

    # make HTTP GET request to path
    headers = {"api-key":"5309619565f744f9248320a886c59bec"}
    r = requests.get(path, params=params, headers=headers)
    
    while r.status_code==503 or r.status_code==502: # Server Error; try again
        attempt += 1
        print(f"Error 503 for {path}; attempt {attempt}", flush=True)
        r = requests.get(path, params=params, headers=headers)

    # raise exception for other response codes that aren't HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = fpath + r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

def periodic_centering(x, center, boxsize):
    quarter = boxsize/4
    upper_qrt = boxsize-quarter
    lower_qrt = quarter
    
    if center > upper_qrt:
        # some of our particles may have wrapped around to the left half 
        x[x < lower_qrt] += boxsize
    elif center < lower_qrt:
        # some of our particles may have wrapped around to the right half
        x[x > upper_qrt] -= boxsize
    
    return x - center

#
# Prep bins, subhalo list
#

# Set bins in r/r200 space
nbins = 128
r_edges = np.logspace(-1, 0.6, nbins+1)
binned_r = r_edges[:-1] + np.diff(r_edges)

folder = "z01_TNG50/"
local = "/mnt/research/galaxies-REU/sims/IllustrisTNG/TNG50/"
littleh = 0.6774
z = 0.01
snapnum = 98
url_dset = "http://www.illustris-project.org/api/TNG50-1/"
url_sbhalos = url_dset + "snapshots/" + str(snapnum) + "/subhalos/"

if rank==0:
    # Get the halos to loop over. It is now "all" of them.
    min_mass = 30* littleh # ~1e11.5 Msun in 1/1e10 Msun / h
    max_mass = 300 * littleh # ~1e12.5 Msun
    search_query = "?mass_dm__gt=" + str(min_mass) \
                 + "&mass_dm__lt=" + str(max_mass)# \
                 #+ "&halfmassrad_stars__gt=" + str(2 / a0 * littleh) # 2 kpc

    cut1 = get(url_sbhalos + search_query)
    cut1['count']
    # Re-get list to avoid paging
    cut1 = get(url_sbhalos + search_query, {'limit':cut1['count'], 'order_by':'id'})

    sub_list = cut1['results']
    sub_ids = np.array([sub['id'] for sub in cut1['results']], dtype='i')

    # Are these subhalos centrals or satellites?
    cat = readsubfHDF5.subfind_catalog(local, snapnum, #grpcat=False, subcat=False,
                                        keysel=['GroupFirstSub','SubhaloGrNr'])
    sat = np.zeros(cat.SubhaloGrNr.size, dtype=bool)
    sat[sub_ids] = (sub_ids != cat.GroupFirstSub[cat.SubhaloGrNr[sub_ids]])
    del cat
    gc.collect()

else:
    sub_ids = None
    sat = None
                                   
my_subs = scatter_work(sub_ids, rank, size)
sub_ids = comm.bcast(sub_ids, root=0)
sat = comm.bcast(sat, root=0)

boxsize = get(url_dset)['boxsize']
a0 = 1/(1+z)

H0 = littleh * 100 * u.km/u.s/u.Mpc

good_ids = np.where(my_subs > -1)[0]
my_profiles = {}

for sub_id in my_subs[good_ids]:

    my_profiles[sub_id] = {}

    #
    # Pull API quantities
    #
    
    sub = get(url_sbhalos + str(sub_id))
    dm_halo = sub["mass_dm"] * 1e10 / littleh * u.Msun
    star_mass = sub["mass_stars"] * 1e10 / littleh * u.Msun
    sfr = sub["sfr"] * u.Msun / u.yr
    r_half = sub["halfmassrad_stars"] * u.kpc * a0/littleh
    x0 = sub["pos_x"]
    y0 = sub["pos_y"]
    z0 = sub["pos_z"]
    
    my_profiles[sub_id]['mass_dm'] = dm_halo
    my_profiles[sub_id]['mass_star'] = star_mass
    my_profiles[sub_id]['ssfr'] = sfr/star_mass
    my_profiles[sub_id]['sat'] = sat[sub_id]

    #
    # Load particle data
    #
    
    gas = True
    
    readhaloHDF5.reset()
    
    try:
        # Gas
        coords = readhaloHDF5.readhalo(local, "snap", snapnum, 
                                        "POS ", 0, -1, sub_id, long_ids=True,
                                        double_output=False).astype("float32")
        
        vel = readhaloHDF5.readhalo(local, "snap", snapnum,
                                    "VEL ", 0, -1, sub_id, long_ids=True,
                                    double_output=False).astype("float32")
        
        dens = readhaloHDF5.readhalo(local, "snap", snapnum, 
                                        "RHO ", 0, -1, sub_id, long_ids=True,
                                        double_output=False).astype("float32")
        
        mass = readhaloHDF5.readhalo(local, "snap", snapnum, 
                                        "MASS", 0, -1, sub_id, long_ids=True,
                                        double_output=False).astype("float32")
        
        inte = readhaloHDF5.readhalo(local, "snap", snapnum, 
                                        "U   ", 0, -1, sub_id, long_ids=True,
                                        double_output=False).astype("float32")
        
        elec = readhaloHDF5.readhalo(local, "snap", snapnum,
                                        "NE  ", 0, -1, sub_id, long_ids=True,
                                        double_output=False).astype("float32")

        cool_rate = readhaloHDF5.readhalo(local, "snap", snapnum,
                                            "GCOL", 0, -1, sub_id, long_ids=True,
                                            double_output=False).astype("float32")
        
    except AttributeError:
        gas = False

    #
    # Calculate r200 and other virial quantities
    #

    r200 = (G*dm_halo/(100*H0**2))**(1/3)
    r200 = r200.to('kpc')

    disp_200 = G*dm_halo/(2*r200)
    disp_200 = np.sqrt(disp_200).to('km/s')
        
    T200 = 0.6*m_p/(2*k_B) * disp_200**2
    T200 = T200.to('K')

    # save virial quantities
    my_profiles[sub_id]['r_200'] = r200
    my_profiles[sub_id]['T_200'] = T200
  
    #
    # Calculate gas information
    #
    
    if gas:
        
        x = coords[:,0]
        y = coords[:,1]
        z = coords[:,2]
        x_rel = periodic_centering(x, x0, boxsize) * u.kpc * a0/littleh
        y_rel = periodic_centering(y, y0, boxsize) * u.kpc * a0/littleh
        z_rel = periodic_centering(z, z0, boxsize) * u.kpc * a0/littleh
        r = np.sqrt(x_rel**2 + y_rel**2 + z_rel**2)
        r_vec = np.column_stack([x_rel,y_rel,z_rel])
        
        mass = mass * 1e10 / littleh * u.Msun
        vel = vel * np.sqrt(a0) * u.km/u.s

        # subtracting off bulk velocity
        vel[:,0] -= sub['vel_x'] * u.km/u.s
        vel[:,1] -= sub['vel_y'] * u.km/u.s
        vel[:,2] -= sub['vel_z'] * u.km/u.s

        #
        # Calculate Entropy
        #

        # For conversion of internal energy to temperature, see
        # https://www.tng-project.org/data/docs/faq/#gen4
        inte *= u.erg/u.g
        X_H = 0.76
        gamma = 5./3.
        mu = 4/(1 + 3*X_H + 4*X_H*elec)
        temp = ( (gamma-1) * inte/k_B * mu*m_p * 1e10 ).to('K')

        dens = dens * 1e10*u.Msun/littleh * (u.kpc*a0/littleh)**-3
        ne = elec * X_H*dens/m_p # elec frac defined as n_e/n_H
        ent = k_B * temp/ne**(gamma-1)
        ent = ent.to('eV cm^2', equivalencies=u.temperature_energy())

        pres = dens/m_p * k_B * temp

        #
        # Calculate some CGM gas properties
        #

        # I will probably change how I distinguish the CGM
        # so that lower radii can be probed
        r_CGM = 2*r_half # DeFelippis+20
        CGM = r > r_CGM
        M_gas_CGM = np.sum(mass[CGM])

        my_profiles[sub_id]['mass_CGM'] = M_gas_CGM

        # Temperatures
        # if ( temp > 1e5*u.K ).any():
        #     hot_CGM = np.logical_and( temp > 1e5*u.K, CGM )
        #     mass_CGM_hot  = np.sum(mass[hot_CGM])
        #     T_hot_avg = np.average(temp[hot_CGM], weights = mass[hot_CGM])
        # else:
        #     mass_CGM_hot  = np.nan * u.Msun
        #     T_hot_avg = np.nan * u.K

        # if ( temp < 1e5*u.K).any():
        #     cool_CGM = np.logical_and( temp < 1e5*u.K, CGM )
        #     mass_CGM_cool = np.sum(mass[cool_CGM])
        # else:
        #     mass_CGM_cool = np.nan * u.Msun
        
        # my_profiles[sub_id]['mass_CGM_hot']  = mass_CGM_hot
        # my_profiles[sub_id]['mass_CGM_cool'] = mass_CGM_cool
        # my_profiles[sub_id]['T_hot_avg'] = T_hot_avg

        #
        # Calculate & store radial profiles
        #
        
        # bin in scaled radial bins
        r_scale = (r/r200).value
        rbinner = np.digitize(r_scale, r_edges)

        binned_ent_avg = np.ones_like(binned_r)*np.nan * u.eV*u.cm**2
        binned_ent_med = np.ones_like(binned_r)*np.nan * u.eV*u.cm**2

        binned_pres_avg = np.ones_like(binned_r)*np.nan * u.dyn/u.cm**2
        binned_pres_med = np.ones_like(binned_r)*np.nan * u.dyn/u.cm**2

        binned_temp_avg = np.ones_like(binned_r)*np.nan * u.K
        binned_temp_med = np.ones_like(binned_r)*np.nan * u.K

        binned_dens_avg = np.zeros_like(binned_r)*np.nan * u.g/u.cm**3
        binned_dens_med = np.zeros_like(binned_r)*np.nan * u.g/u.cm**3
        
        # find central tendency for each radial bin
        for i in range(1, r_edges.size):
            this_bin = rbinner==i
            if np.sum(mass[this_bin]) != 0: # are there particles in this bin

                binned_ent_avg[i-1] = np.average(ent[this_bin],
                                                 weights = mass[this_bin])
                binned_ent_med[i-1] = np.median(ent[this_bin])

                binned_pres_avg[i-1] = np.average(pres[this_bin],
                                                  weights = mass[this_bin])
                binned_pres_med[i-1] = np.median(pres[this_bin])

                binned_temp_avg[i-1] = np.average(temp[this_bin],
                                                  weights = mass[this_bin])
                binned_temp_med[i-1] = np.median(temp[this_bin])

                binned_dens_avg[i-1] = np.average(dens[this_bin],
                                                  weights = mass[this_bin])
                binned_dens_med[i-1] = np.median(dens[this_bin])
                
        my_profiles[sub_id]['ent_avg'] = binned_ent_avg
        my_profiles[sub_id]['ent_med'] = binned_ent_med
        my_profiles[sub_id]['pres_avg'] = binned_pres_avg
        my_profiles[sub_id]['pres_med'] = binned_pres_med
        my_profiles[sub_id]['temp_avg'] = binned_temp_avg
        my_profiles[sub_id]['temp_med'] = binned_temp_med
        my_profiles[sub_id]['dens_avg'] = binned_dens_avg
        my_profiles[sub_id]['dens_med'] = binned_dens_med

    else: # no gas
        my_profiles[sub_id]['mass_CGM'] = np.nan * u.Msun
        # my_profiles[sub_id]['mass_CGM_hot']  = np.nan * u.Msun
        # my_profiles[sub_id]['mass_CGM_cool'] = np.nan * u.Msun
        # my_profiles[sub_id]['T_hot_avg'] = np.nan * u.K
        
        my_profiles[sub_id]['ent_avg'] = np.nan
        my_profiles[sub_id]['ent_med'] = np.nan
        my_profiles[sub_id]['pres_avg'] = np.nan
        my_profiles[sub_id]['pres_med'] = np.nan
        my_profiles[sub_id]['temp_avg'] = np.nan
        my_profiles[sub_id]['temp_med'] = np.nan
        my_profiles[sub_id]['dens_avg'] = np.nan
        my_profiles[sub_id]['dens_med'] = np.nan
        
#
# Collect data from MPI ranks & write to files
#

profile_list = comm.gather(my_profiles, root=0)

if rank==0:

    all_gal_prop = np.zeros( (len(sub_ids), 9) )
    all_ent_prof = np.zeros( (len(sub_ids), 2*nbins+1) )
    all_pres_prof = np.zeros( (len(sub_ids), 2*nbins+1) )
    all_temp_prof = np.zeros( (len(sub_ids), 2*nbins+1) )
    all_dens_prof = np.zeros( (len(sub_ids), 2*nbins+1) )
    
    i=0
    for dic in profile_list:
        for k,v in dic.items():
            all_gal_prop[i,0] = k
            all_gal_prop[i,1] = v['sat']
            all_gal_prop[i,2] = v['mass_dm'].value
            all_gal_prop[i,3] = v['mass_star'].value
            all_gal_prop[i,4] = v['mass_CGM'].value
            all_gal_prop[i,5] = v['ssfr'].value
            all_gal_prop[i,6] = v['r_200'].value
            all_gal_prop[i,8] = v['T_200'].value
            # all_gal_prop[i,9] = v['mass_CGM_hot'].value
            # all_gal_prop[i,10] = v['mass_CGM_cool'].value
            # all_gal_prop[i,11] = v['T_hot_avg'].value

            all_ent_prof[i,0] = k
            all_ent_prof[i,1::2] = v['ent_avg']
            all_ent_prof[i,2::2] = v['ent_med']
            
            all_pres_prof[i,0] = k
            all_pres_prof[i,1::2] = v['pres_avg']
            all_pres_prof[i,2::2] = v['pres_med']
            
            all_temp_prof[i,0] = k
            all_temp_prof[i,1::2] = v['temp_avg']
            all_temp_prof[i,2::2] = v['temp_med']

            all_dens_prof[i,0] = k
            all_dens_prof[i,1::2] = v['dens_avg']
            all_dens_prof[i,2::2] = v['dens_med']

            i+=1

    sort = np.argsort(all_gal_prop[:,0])

    prop_header = "SubID,Sat,MassDark,MassStar,MassCGM,sSFR,r200,T200,THot"

    header = "SubID"
    for r in binned_r:
        header += "   {:.4f} avg med".format(r)

    np.savetxt(folder+'halo_properties.csv', all_gal_prop[sort],
               delimiter=',', header=prop_header)
    np.savetxt(folder+'profiles_entropy.csv', all_ent_prof[sort], 
               delimiter=',', header=header)
    np.savetxt(folder+'profiles_pressure.csv', all_pres_prof[sort],
               delimiter=',', header=header)
    np.savetxt(folder+'profiles_temperature.csv', all_temp_prof[sort], 
               delimiter=',', header=header)
    np.savetxt(folder+'profiles_density.csv', all_dens_prof[sort],
               delimiter=',', header=header)
