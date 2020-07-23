import numpy as np
import os
from scipy.special import \
	hyp2f1
from pygrackle import \
    chemistry_data, \
    setup_fluid_container
from pygrackle.utilities.physical_constants import \
	cm_per_mpc as Mpc, \
	amu_cgs as m_u, \
	mass_sun_cgs as Msun, \
	keV_per_K as kboltz, \
	sec_per_Myr, \
	cm_per_kpc
	
global Gamma
global mu

Gamma     = 5.0/3.0
mu        = 0.6

def mass_integral(R_vir, r_s, a, S_0, T_kev, S_c=0):
	#S_0 = T_kev/n_0**(Gamma-1)
	
	if S_c != 0:
		mass = ( S_c + S_0*(R_vir/r_s)**a )**(1/(1-Gamma)) \
		     * ( T_kev*(1 + S_0*(R_vir/r_s)**a / S_c) )**(1/(Gamma-1)) \
		     * hyp2f1(3/a, 1/(Gamma-1), 1+3/a, -S_0*(R_vir/r_s)**a/S_c) \
		     * 4.0/3.0*np.pi*mu*m_u * R_vir**3
		
	else:
		mass = -4*np.pi*mu*m_u * R_vir**3 * (S_0*(R_vir/r_s)**a)**(1/(1-Gamma))\
		     * (Gamma - 1) * T_kev**(1/(Gamma-1))
		mass /= 3+a-3*Gamma
		
	mass /= Msun
	return mass

def radial_bins(R_vir, r_s, a, n_0, S_0, T_kev, S_c=0, nbins=8192):
	r = np.linspace(0, R_vir, nbins)
	dr = r[2]-r[1]
	temp = T_kev/kboltz*np.ones(nbins) # K
	
	if S_c != 0:
		entr = S_c + S_0*(r/r_s)**a
		dens = mu*m_u * (T_kev / (S_c + S_0*(r/r_s)**a))**(1/(Gamma-1)) # g/cm^3
		dens[0] = 0
		
	else:
		entr = S_0*(r/r_s)**a
		dens = mu*m_u * n_0*(r/r_s)**(-a/(Gamma-1)) # g/cm^3
		dens[0] = 0
	
	mass = 4/3*np.pi*dens*((r+dr)**3 - r**3) # cgs
	mass /= Msun
	
	return r, temp, entr, dens, mass

# Parameter file values
GalaxySimulationGasHaloTemperature = 1.5e6
GalaxySimulationGasHaloDensity = 1.5e-27
GalaxySimulationGasHaloScaleRadius = 0.06
GalaxySimulationGasHaloAlpha = 2.0/3.0
GalaxySimulationGasHaloCoreEntropy = 0 # must be set

M = 1e12 * Msun
rho_crit = 1.8788e-29*0.49
R_vir = (3.0/(4.0*3.14159)*M/(200.*rho_crit)) ** (1./3.)

# Devin's values
D_a = GalaxySimulationGasHaloAlpha # rename
D_r_s = GalaxySimulationGasHaloScaleRadius * Mpc # scale radius in cm
D_T_kev = GalaxySimulationGasHaloTemperature * kboltz  # halo temperature in keV
D_n_0 = GalaxySimulationGasHaloDensity / (mu*m_u)  # convert n_0 to electron number density
D_S_0 = D_T_kev / pow(D_n_0,Gamma-1.0)
D_S_c = GalaxySimulationGasHaloCoreEntropy # rename

#print(D_S_0)

# Parameter adjustments
S_0 = D_S_0 # keep this
S_c = D_S_c # and this
a = D_a # and this

r_s = D_r_s#0.05 * Mpc # keep this, at least for now

#n_0 = D_n_0
n_0 = 0.5e-27 #0.9e-27 # set as density in cgs
n_0 /= mu * m_u # to number density

#T_kev = D_T_kev
T_kev = S_0 * pow(n_0,Gamma-1.0)
T_K = T_kev / kboltz

mass_int = mass_integral(R_vir, r_s, a, S_0, T_kev, S_c=S_c) 

r, temp, entr, dens, mass = radial_bins(R_vir, r_s, a, n_0, S_0, T_kev, S_c=S_c)

print("{:.3e} Msun".format(mass_int))
print("{:.3e} Msun".format(mass.sum()))
print('{:.3e} K'.format(T_K))

current_redshift = 0.

# Set solver parameters
chem = chemistry_data()
chem.use_grackle = 1
chem.with_radiative_cooling = 1
chem.primordial_chemistry = 2
chem.metal_cooling = 1
chem.UVbackground = 1
chem.cmb_temperature_floor = 1
chem.h2_on_dust = 1
chem.grackle_data_file = "grackle/input/CloudyData_UVB=HM2012.h5"

chem.use_specific_heating_rate = 1
chem.use_volumetric_heating_rate = 1

# Set units
chem.comoving_coordinates = 0 # proper units
chem.a_units = 1.0
chem.a_value = 1.0 / (1.0 + current_redshift) / chem.a_units
chem.density_units = m_u # rho = 1.0 is ~1.67e-24 g/cm^3
chem.length_units = 1638.4 * cm_per_kpc
chem.time_units = sec_per_Myr
chem.velocity_units = chem.a_units \
    * (chem.length_units / chem.a_value) \
    / chem.time_units

# Call convenience function for setting up a fluid container.
# This container holds the solver parameters, units, and fields.
fc = setup_fluid_container(chem,
                           density=dens[1:].copy(),
                           temperature=temp[0],#.copy(),
                           metal_mass_fraction=0.002041,
                           converge=False)

#fc["specific_heating_rate"][:] = 0.
#fc["volumetric_heating_rate"][:] = 0.

#fc.calculate_temperature() # update to better reflect energy
fc.calculate_cooling_time() # Myr; negative means cooling
print(fc['cooling_time'][(r[:-1]/cm_per_kpc) > 25]/sec_per_Myr)
