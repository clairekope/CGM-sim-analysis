import numpy as np
import unyt as u
import unyt.dimensions as udim
from unyt import accepts, returns

class DiskModel():

    @accepts(Mvir=udim.mass, Mstar=udim.mass, rs_star=udim.length, zs_star=udim.length)
    def __init__(self, 
                 M_vir=1e12*u.Msun,
                 C=10, 
                 M_star=5.8e10*u.Msun,
                 rs_star=3.5*u.kpc,
                 zs_star=0.325*u.kpc):
        """
        Class for calculating the Toomre criterion as a function of radius.
        Initializer sets some basic parameters about the disk model.

        Mvir (unyt quantity)
            Virial (dark matter) mass of the system.
            Default: 1e12 Msun

        C (int or float)
            Concentration of the NFW profile.
            Default: 10

        Mstar (unyt quantity)
            Stellar mass for the disk potential.
            Default: 5.8e10 Msun

        rs_star (unyt quantity)
            Scale radius for the stellar disk potential.
            Default: 3.5 kpc

        zs_star (unyt quantity)
            Scale height for the stellar disk potential.
            Default: 0.325 kpc
        """
        
        self.M_vir = M_vir
        self.C = C
        self.M_star = M_star
        self.rs_star = rs_star
        self.zs_star = zs_star

        self.rho_crit = 1.8788e-29*0.49 * u.g/u.cm**3
        self.Rvir = ( 3.0/(4.0*np.pi)*self.M_vir / (200.*self.rho_crit) )**(1./3.)
        self.Rs = self.Rvir/self.C
        self.M_NFW = self.M_vir /  (np.log(1.0+self.C) - self.C/(1.0+self.C))

    @accepts(radii=udim.length)
    @returns(udim.length/udim.time**2)
    def g(self, radii):
        """
        Gravitational acceleration in the cylindrical r direction for z=0 (aka, in the disk plane).
        A combination of NFW and Miyamoto & Nagai potentials. Ignores a bulge potential.

        Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration

        Outputs:
        --------
        g_r (unyt array)
            Radial accelerations at each input radius
        """

        # NFW potential
        x = radii/self.Rs

        g_NFW = u.G*self.M_NFW/radii**2 * (np.log(1+x) - x/(1+x))
        
        # Miyamoto & Nagai potential
        k = self.rs_star + self.zs_star

        g_MN = u.G*self.M_star*radii / np.power( np.power(radii,2) + k**2, 3/2 )

        return g_NFW + g_MN

    @accepts(radii=udim.length)
    @returns(1/udim.time**2)
    def dg_dr(self, radii):
        """
        Radial derivative of the gravitational acceleration in the cylindrical r direction, at z=0 (aka, in the disk plane).
        A combination of NFW and Miyamoto & Nagai potentials. Ignores a bulge potential.

        Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration's derivative

        Outputs:
        --------
        dg_dr (unyt array)
            Derivative of the radial acceleration at each input radius
        """

        # NFW potential
        x = radii/self.Rs

        dg_NFW = np.power(radii, -2) * (2*self.Rs +3*radii)/np.power(self.Rs+radii,2)
        dg_NFW -= 2*np.power(radii,-3)*np.log(1+x)
        dg_NFW *= u.G * self.M_NFW

        # Miyamoto & Nagai potential
        k = self.rs_star + self.zs_star
        d = np.power(radii, 2) + k**2

        dg_MN = u.G*self.M_star * np.power(d, -3/2) * (1 - 3*np.power(radii, 2)/d)

        return dg_NFW + dg_MN

    @accepts(radii=udim.length)
    @returns(1/udim.time)
    def omega(self, radii):
        """
        Rotation frequency Omega = v_phi/r, where v_phi = sqrt(r * g(r))

        Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration's derivative

        Outputs:
        --------
        omega (unyt array)
            Rotation frequency as a function of radius. Radians are implied.
        """

        return np.sqrt(self.g(radii) / radii)

    @accepts(radii=udim.length)
    @returns(udim.length/udim.time)
    def _epicyclic_derivative(self, radii):
        """
        Calculates the d/dr (r^2 * Omega) term of the epicyclic frequency, kappa.

        Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration's derivative

        Outputs:
        --------
        ep_deriv (unyt array)
            Derivative of r^2 * Omega at supplied radii
        """

        g = self.g(radii)
        dg_dr = self.dg_dr(radii)

        num = 3*np.power(radii,2)*g + np.power(radii,3)*dg_dr
        denom = 2*np.sqrt( np.power(radii, 3) * g )

        return num/denom

    @accepts(radii=udim.length)
    @returns(1/udim.time)
    def kappa(self, radii):
        """
        Epicyclic frequency as a function of radius

        Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration's derivative

        Outputs:
        --------
        kappa (unyt array)
            Epicyclic frequency at supplied radii
        """

        omega = self.omega(radii)
        deriv = self._epicyclic_derivative(radii)

        return np.sqrt(2*omega/radii * deriv)

if __name__ == "__main__":

    r = np.linspace(0,30) * u.kpc
    dsk = DiskModel()

    # Function decorators ensure return products have desired units
    dsk.g(r)
    dsk.dg_dr(r)
    dsk.omega(r)
    dsk._epicyclic_derivative(r)
    dsk.kappa(r)