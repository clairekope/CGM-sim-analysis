# "translated" and abridged from ESDescriptor-Greg Meece

from __future__ import division, print_function
import numpy as np
import pdb
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings('error')

#-----Setup-----
#r_min = 1.0e-3 # Mpc
#r_max = 1.0 # Mpc
r_min = 1e-6
r_max = 2000.0
entropy_k0 = 4.0 # keV cm^2
entropy_k1 = 20.0 # keV cm^2
entropy_slope = 1.15
entropy_scale_length = 0.1 # Mpc
nfw_mvir = 6.0e+11 # Param file; 1e12?
nfw_concentration = 12.0 # Param file; GalaxySimulationDarkMatterConcentrationParameter

length_units = 3.08567758e+21 # 1 kpc in cm
dens_units = 1.66054e-24  # 1 particle/cc in g cm^-3
time_units = 3.15576e+14 # 10 Million years in s
mass_units = dens_units * pow(length_units, 3.0);
velocity_units = length_units / time_units;

mu = 0.6 # ionized gas
G_CGS = 6.67259e-8
AMU_CGS = 1.6605402e-24
KB_CGS = 1.380658e-16
CM_PER_KM = 1.0e5
CM_PER_MEGAPARSEC = 3.085677581e24
VIRIAL_COEFFICIENT = 200.0
SOLAR_MASS_IN_GRAMS = 1.9891e+33
ERG_PER_KEV = 1.6021773e-9
Gamma = 5.0/3.0
H0 = 70.0
omega_m = 0.3
omega_l = 0.7
redshift = 0.0
gravitational_constant_cgs = 6.67259e-8

mp_code = AMU_CGS / mass_units
kb_code = KB_CGS / (mass_units * velocity_units * velocity_units)
g_code = gravitational_constant_cgs * dens_units * pow(time_units, 2.0)

hubble_constant = H0 * np.sqrt(omega_m * pow(1.0 + redshift, 3.0) + omega_l) # km/s/mpc
hubble_constant *= time_units * CM_PER_KM / CM_PER_MEGAPARSEC # code
critical_density = 3.0 * pow(hubble_constant, 2.0) / (8.0 * np.pi * g_code) # code
nfw_mvir *= SOLAR_MASS_IN_GRAMS / (dens_units * pow(length_units, 3.0)) # code
nfw_scale_radius = 3.0*nfw_mvir / (4.0*np.pi*VIRIAL_COEFFICIENT*critical_density*pow(nfw_concentration, 3.0)) 
nfw_scale_radius = pow(nfw_scale_radius, 1.0/3.0)
nfw_mass_temperature = (mu * mp_code/(2.0*kb_code)) * pow(10.0*g_code*nfw_mvir*hubble_constant, 2.0/3.0) #code 

print("Scale radius:", nfw_scale_radius)

#-----Functions-----
# Integrates bin by bin
def integrate_density(r0,r1,rho):
	rho_orig = rho
	max_dr = 1.0e-4
	#Integrating outwards
	if r0 < r1:
		print(r0)
		r = r0
		while r < r1:
			dr = min(max_dr, r1 - r)				
			rho = rk4(r, rho, dr)
			r += dr

	#Integrating inwards
	if r0 > r1:
		r = r0
		while r > r1:
			dr = min(max_dr, r - r1)
			dr = -dr			
			# This isn't yet causing problems so we'll leave it
			rho = rk4(r, rho, dr)
			r += dr
	return rho

def rk4(r,rho,dr):
	k1 = drhodz(r, rho)
	k2 = drhodz(r + 0.5*dr, rho + 0.5*k1*dr)
	k3 = drhodz(r + 0.5*dr, rho + 0.5*k2*dr)
	k4 = drhodz(r + dr, rho + k3*dr)

	return rho + (dr/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)

def drhodz(r,rho):
	g = g_code * get_total_mass(r) / pow(r, 2.0); # code
	g *= length_units / pow(time_units, 2.0); # cgs
	
	r *= length_units # cgs
	rho *= dens_units # cgs
	esl = entropy_scale_length * length_units # cgs
	
	k0 = entropy_k0 * ERG_PER_KEV # cgs
	k1 = entropy_k1 * ERG_PER_KEV # cgs
	k = k0 + k1 * pow(r/esl, entropy_slope) # cgs
	dkdr = entropy_slope * (k1 /esl) * pow(r/esl, entropy_slope - 1.0) # cgs

	fh = 0.76
	alpha = (fh + 2.0*(1.0-fh)/4.0) * mu  # ne = alpha * n
	n = rho / (mu * AMU_CGS) # number density
	#ne = alpha * n # Electron number density

	drho = -( mu * AMU_CGS * n * g / pow(alpha, Gamma - 1.0) + dkdr * pow(n, Gamma) )
	drho /= Gamma * pow(n, Gamma - 1.0) * k # dndr in cgs
	drho *= mu * AMU_CGS # drhodr in cgs
	drho /= (dens_units / length_units) # code

	return drho

def get_total_mass(r): # get_nfw_mass 
	# Take r in code units
	m_tot = np.log( (r+nfw_scale_radius)/nfw_scale_radius ) - r/(r+nfw_scale_radius)
	m_tot /= np.log(1.0 + nfw_concentration) - nfw_concentration/(1.0+nfw_concentration)
	m_tot *= nfw_mvir # code

	return m_tot

#-----Initialize bins-----
#Greg use 512 radial bins
n_radial_bins = 512
bin_left = np.zeros(n_radial_bins)
bin_dens = np.zeros(n_radial_bins)
bin_temp = np.zeros(n_radial_bins)

#-----Setup radial bins----- 
#r_min *= CM_PER_MEGAPARSEC # cgs
#r_max *= CM_PER_MEGAPARSEC # cgs
#r_min /= length_units # code (kpc, in this case)
#r_max /= length_units # code

log_r_min = np.log10(r_min)
log_r_max = np.log10(r_max)

for i in range(n_radial_bins):
	bin_left[i] = log_r_min + (i/n_radial_bins) * (log_r_max-log_r_min)
	bin_left[i] = pow(10.0, bin_left[i])

#Find radial bin that holds rv the viral radius
rv = nfw_scale_radius * nfw_concentration
#print(rv)

if rv < bin_left[0]: #less than r_min
	vbin = -1
elif rv > bin_left[n_radial_bins - 1]: #greater than r_max
	vbin = n_radial_bins - 1
else:
	for i in range(n_radial_bins-1):
		if bin_left[i] <= rv and bin_left[i+1] > rv:
			vbin = i

#-----Find Density-----
entropy_scale_length *= CM_PER_MEGAPARSEC # cgs
entropy_scale_length /= length_units # code

#Calculate density of virial bin/radius from entropy
#Compute in cgs first, then convert 
vk = entropy_k0 + entropy_k1 * pow(rv/entropy_scale_length, entropy_slope) # KeV cm^2
vk *= ERG_PER_KEV # erg cm^2

fh = 0.76;
vdens = KB_CGS * nfw_mass_temperature / vk # ne^(2/3)
vdens = pow(vdens, 3.0/2.0) # ne
vdens /= mu * (fh + (1.0-fh)/2.0) # n
vdens *= mu * AMU_CGS # rho in cgs
vdens /= dens_units #rho in code

#Set vbin, accounting that rv is not a bin edge
if vbin == -1: # rv < r_min
	bin_dens[0] = integrate_density(rv, bin_left[0], vdens);
	vbin = 0
else:
	bin_dens[vbin] = integrate_density(rv, bin_left[vbin], vdens)

#Integrate to find the density in other bins
	for i in range(vbin-1,-1,-1):
		bin_dens[i] = integrate_density(bin_left[i+1], bin_left[i], bin_dens[i+1])

	for i in range(vbin,n_radial_bins - 1):
		bin_dens[i+1] = integrate_density(bin_left[i], bin_left[i+1], bin_dens[i])

#print(bin_dens)

#-----Find Temperature-----
for i in range(n_radial_bins):
	ki = entropy_k0 + entropy_k1 * pow(bin_left[i] / entropy_scale_length, entropy_slope) # Kev cm^2
	ki *= ERG_PER_KEV # Ergs cm^2
	n = (bin_dens[i] * dens_units) / (mu * AMU_CGS)
	fh = 0.76
	n *= (fh + (1.0-fh)/2.0) * mu #ne

	bin_temp[i] = ki * pow(n, 2.0/3.0) / KB_CGS
	
#-----Print & Plot-----
#for i in range(n_radial_bins):
#  print("{3d} {10.8f} {10.8f} {10.8f} {10.8f}".format(i, bin_left[i], bin_temp[i], bin_dens[i]))

plt.loglog(bin_left,bin_dens)
plt.loglog(bin_left,bin_temp)
plt.legend(["Density","Temperature"],loc='lower left')
plt.show()
 
