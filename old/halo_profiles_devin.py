import numpy as np

GalaxySimulationGasHaloGalaxyMass      = 1e12
GalaxySimulationGasHaloTemperature     = 1.5e6
GalaxySimulationGasHaloDensity         = 1.5e-27
GalaxySimulationGasHaloCoreEntropy     = 5.0
GalaxySimulationGasHaloAlpha           = 0.6667
GalaxySimulationGasHaloScaleRadius     = 0.06 # For entropy power law
GalaxySimulationGasHaloDMConcentration = 12

SolarMass = 1.99e33 # cgs
GravConst = 6.67e-8 
Mpc	      = 3.0856e24 
Gamma     = 5/3 # pretty sure
kboltz    = 1.381e-16
rho_crit  = 1.8788e-29*0.49
mh        = 1.67e-24

class CGMData():
    def __init__(self, nbins, R_outer):
    
        self.nbins = nbins
        
        self.n_rad = -1*np.ones(self.nbins)
        self.T_rad = -1*np.ones(self.nbins)
        self.rad   = -1*np.ones(self.nbins)

        self.R_outer = R_outer
        self.dr = self.R_outer/self.nbins

def halo_S_of_r(r):
    """
    Halo entropy as a function of radius for the user-specified CGM types that require numerical
    integration. 

    Input is radius in CGS units.  output is entropy in CGS units (Kelvin cm^2)
    """

    # calculate a bunch of things based on user inputs
    Tvir = GalaxySimulationGasHaloTemperature  # in Kelvin
    n0 = GalaxySimulationGasHaloDensity / (0.6*1.67e-24)  # convert from density to electron number density (cm^-3)
    r0 = GalaxySimulationGasHaloScaleRadius * Mpc  # scale radius in CGS
    Smin = GalaxySimulationGasHaloCoreEntropy/8.621738e-8  # given in keV cm^2, converted to Kelvin cm^2
    S0 = Tvir / pow(n0,Gamma-1)  # entropy at scale radius, in units of Kelvin cm^2

    return Smin + S0*pow(r/r0,GalaxySimulationGasHaloAlpha)  # has units of Kelvin cm^2
    #return S0*pow(r/r0,GalaxySimulationGasHaloAlpha)

def halo_dSdr(r):

	"""
	dEntropy/dr as a function of radius for the user-specified CGM types that require numerical
	integration. 

	Input is radius in CGS units output is entropy gradient in CGS units (Kelvin*cm)
	"""

	# calculate a bunch of things based on user inputs
	Tvir = GalaxySimulationGasHaloTemperature  # in Kelvin
	n0 = GalaxySimulationGasHaloDensity / (0.6*1.67e-24)  # convert from density to electron number density (cm^-3)
	r0 = GalaxySimulationGasHaloScaleRadius*Mpc  # scale radius in CGS
	Smin = GalaxySimulationGasHaloCoreEntropy/8.621738e-8  # given in keV cm^2, converted to Kelvin cm^2
	S0 = Tvir / pow(n0,Gamma-1)  # entropy at scale radius, in units of Kelvin cm^2

	# has units of Kelvin*cm, same deriv for both halo types (since constant drops out)
	return S0*GalaxySimulationGasHaloAlpha * \
	pow(r/r0,GalaxySimulationGasHaloAlpha-1.0)/r0

def halo_dn_dr(r, n):
    """
    dn/dr as a function of radius and halo electron number density.  This quantity is calculated
    by assuming that gravity and pressure are in hydrostatic equilibrium in a halo with a specified 
    entropy profile S(r).

    Input is radius in CGM units and electron number density in units 
    of particles per cm^-3.  Output is dn/dr in CGS units, so particles per cm^4. 
    """

    return -1.0*( n*1.22*mh*halo_g_of_r(r) + kboltz*pow(n,Gamma)*halo_dSdr(r) ) \
    / ( Gamma * kboltz * halo_S_of_r(r) * pow(n, Gamma-1))

	
def halo_g_of_r(r):
	"""
	halo gravitational acceleration as a function of radius.

	Input is the radius in CGS units and returns the MAGNITUDE of the 
	acceleration in CGS units.  
	"""
	return GravConst*halo_galmass_at_r(r)/(r*r) 

def halo_galmass_at_r(r):
	"""
	halo galaxy mass at a given radius, using user-defined global parameters for galaxy
	quantities and assuming that all halo mass is in an NFW halo.  This is not totally
	correct near the center of the halo, but since we're using it for the CGM initialization 
	and are dealing with radii that aren't particularly near the center of the halo, this 
	approximation is probably fine. 

	Input is the radius in CGS units output is the enclosed mass at that radius in CGS units.
	"""

	M = GalaxySimulationGasHaloGalaxyMass * SolarMass  # halo total mass in CGS
	C = GalaxySimulationGasHaloDMConcentration  # concentration parameter for NFW halo

	Rvir = pow(3.0/(4.0*np.pi)*M/(200.*rho_crit),1./3.)  # virial radius in CGS
	Rs = Rvir/C  # scale radius of NFW halo in CGS
	rho_0 = 200.0*pow(C,3)/3.0/(np.log(1.0+C) - C/(1.0+C))*rho_crit  # rho_0 for NFW halo in CGS

	# mass w/in radius R
	M_within_r = 4.0*np.pi*rho_0*pow(Rs,3.0)*(np.log((Rs+r)/Rs) - r/(Rs+r))

	return M_within_r

########
# MAIN #
########

M = GalaxySimulationGasHaloGalaxyMass * SolarMass  # halo total mass in CGS
Rvir = pow(3.0/(4.0*np.pi)*M/(200.*rho_crit),1./3.)  # virial radius in CGS

print(Rvir)

CGM_data = CGMData(8192, Rvir)

# set some quantities based on user inputs this defines our integration
Tvir = GalaxySimulationGasHaloTemperature
n0 = GalaxySimulationGasHaloDensity / (0.6*1.67e-24)
r0 = GalaxySimulationGasHaloScaleRadius * Mpc

# used for our numerical integration
dr = CGM_data.dr
this_n = n0
this_radius = r0

# set the bin that we start at (otherwise it doesn't get set!)
index = int(this_radius/dr+1.0e-3)
CGM_data.n_rad[index] = this_n
CGM_data.T_rad[index] = Tvir
CGM_data.rad[index] = this_radius

# starting at the point where the user has defined the radius, density, and 
# temperature, use RK4 to integrate the number density outward to R_outer using the expression 
# for dn_dr in another function.  Calculate the temperature using the entropy at this radius. 
while this_radius <= CGM_data.R_outer-dr:
    #print("integrate outward")
    #print(this_radius, CGM_data.R_outer)

    # calculate RK4 coefficients.
    k1 = halo_dn_dr(this_radius, this_n)
    k2 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k1)
    k3 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k2)
    k4 = halo_dn_dr(this_radius + dr, this_n + dr*k3)

    # update density and radius
    this_n += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4)
    this_radius += dr  # new radius

    # calculate temperature at this radius using known entropy 
    temperature = halo_S_of_r(this_radius) * pow(this_n,Gamma-1.0)

    # store everything in the struct
    index = int(this_radius/dr+1.0e-3)    
    CGM_data.n_rad[index] = this_n
    CGM_data.T_rad[index] = temperature
    CGM_data.rad[index] = this_radius

# now we do the same thing as above, but integrating intward to zero radius.
this_n = n0
this_radius = r0
dr *= -1.0

while this_radius > -dr:
    #print("integrate inward")

    # calculate RK4 coefficients.
    k1 = halo_dn_dr(this_radius, this_n)
    k2 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k1)
    k3 = halo_dn_dr(this_radius + 0.5*dr, this_n + 0.5*dr*k2)
    k4 = halo_dn_dr(this_radius + dr, this_n + dr*k3)

    # update density and radius
    this_n += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4)
    this_radius += dr  # new radius

    # calculate temperature at this radius using known entropy 
    temperature = halo_S_of_r(this_radius) * pow(this_n,Gamma-1.0)

    # store everything in the struct
    index = int(this_radius/(-1.0*dr)+1.0e-3)

    if index >= 0:
      CGM_data.n_rad[index] = this_n
      CGM_data.T_rad[index] = temperature
      CGM_data.rad[index] = this_radius

# this integration acts a little squirrelly around r=0 because the mass values are garbage.  Cheap fix.
#CGM_data.rad[0]=CGM_data.rad[1]
#CGM_data.n_rad[0]=CGM_data.n_rad[1]
#CGM_data.T_rad[0]=CGM_data.T_rad[1]
