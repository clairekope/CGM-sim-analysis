import numpy as np
import unyt as u
import unyt.dimensions as udim
from unyt import accepts, returns

class DiskModel():

    @accepts(Mvir=udim.mass, Mstar=udim.mass,
             M_disk=udim.mass, T_disk=udim.temperature,
             r_scale=udim.length, z_scale=udim.length)
    def __init__(self, 
                 M_vir=1e12*u.Msun,
                 C=10, 
                 M_star=5.8e10*u.Msun,
                 M_disk=5e9*u.Msun,
                 T_disk=1e5*u.K,
                 r_scale=3.5*u.kpc,
                 z_scale=0.325*u.kpc,):
        """
        Class for calculating the Toomre parameter as a function of radius.
        Initializer sets some basic parameters about the disk model.

        M_vir (unyt quantity)
            Virial (dark matter) mass of the system.
            Default: 1e12 Msun

        C (int or float)
            Concentration of the NFW profile.
            Default: 10

        M_star (unyt quantity)
            Stellar mass for the disk potential.
            Default: 5.8e10 Msun

        M_disk (unyt quantity)
            Mass of the gas disk used for defining the gas density profile.
            Default: 5e9 Msun

        T_disk (unyt quantity)
            Temperature of the isothermal gaseous disk. Assumed to be an ideal gas.
            Default: 1e5 K

        r_scale (unyt quantity)
            Scale radius for the stellar disk potential & gas density.
            Default: 3.5 kpc

        z_scale (unyt quantity)
            Scale height for the stellar disk potential & gas density.
            Default: 0.325 kpc
        """
        
        self.M_vir = M_vir
        self.C = C
        self.M_star = M_star
        self.M_disk = M_disk
        self.T_disk = T_disk
        self.r_scale = r_scale
        self.z_scale = z_scale

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
        k = self.r_scale + self.z_scale

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

        dg_NFW = np.power(radii, -2) * (2*self.Rs + 3*radii)/np.power(self.Rs+radii,2)
        dg_NFW -= 2*np.power(radii,-3)*np.log(1+x)
        dg_NFW *= u.G * self.M_NFW

        # Miyamoto & Nagai potential
        k = self.r_scale + self.z_scale
        d = np.power(radii, 2) + k**2

        dg_MN = u.G*self.M_star * np.power(d, -3/2) * (1 - 3*np.power(radii, 2)/d)

        return dg_NFW + dg_MN

    @accepts(radii=udim.length)
    @returns(1/udim.time)
    def Omega(self, radii):
        """
        Rotation frequency Omega = v_phi/r, where v_phi = sqrt(r * g(r))

        Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration's derivative

        Outputs:
        --------
        Omega (unyt array)
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

        Omega = self.Omega(radii)
        deriv = self._epicyclic_derivative(radii)

        return np.sqrt(2*Omega/radii * deriv)

    @accepts(radii=udim.length)
    @returns(udim.length/udim.time)
    def c_sound(self, radii):
        """
        Calculate the sound speed of an isothermal disk (i.e., constant sound speed), 
        assuming an ideal gas and mean molecular weight = 0.6 (primordial).

        Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration's derivative

        Outputs:
        --------
        cs (unyt array)
            Sound speed at supplied radii
        """

        gamma = 5/3

        return np.sqrt(gamma * u.kb * self.T_disk / (0.6*u.mp))

    @accepts(radii=udim.length)
    @returns(udim.mass/udim.length**2)
    def Sigma(self, radii):
        """
        Calculate the surface density of the disk profile in Tonnesen & Bryan 09.
        This is a double-sech disk with a smoothing factor applied for cylindrical radii r > 24 kpc.
        The double-sech profile makes it easy to integrate out the vertical (z) component.

         Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration's derivative

        Outputs:
        --------
        Sigma (unyt array)
            Surface density at the supplied radii
        """

        R = 24*u.kpc

        surdens = self.M_disk/(8*self.r_scale**2) / np.cosh(radii/self.r_scale)
        smooth = 0.5 * (1 + np.cos( np.pi * (radii-R)/(7.2*u.kpc) ))

        # multiply radii > R by smoothing factor
        # np.where strips units so reapply them
        return np.where(radii < R, surdens, surdens * smooth) * surdens.units

    @accepts(radii=udim.length)
    @returns(udim.dimensionless)
    def Q(self, radii):
        """
        Calculate the Toomre parameter for an isothermal disk of ideal gas.

         Inputs:
        -------
        radii (unyt array)
            Radii at which to calculate the acceleration's derivative

        Outputs:
        --------
        Q (unyt array)
            The Toomre parameter
        """

        kappa = self.kappa(radii) # Depends on the preceeding functions
        cs = self.c_sound(radii)
        Sigma = self.Sigma(radii)

        return cs*kappa/(np.pi*u.G*Sigma)

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    r = np.linspace(0,30) * u.kpc

    # Function decorators ensure return products have desired units!
    dsk = DiskModel()
    Q_hot = dsk.Q(r)
    
    dsk2 = DiskModel(T_disk=1e4*u.K)
    Q_cool = dsk2.Q(r)

    dsk3 = DiskModel(T_disk=1e3*u.K)
    Q_cold = dsk3.Q(r)

    plt.semilogy(r, Q_hot, label="T=1e5 K", color='C1')
    plt.semilogy(r, Q_cool, label="T=1e4 K", color='C2')
    plt.semilogy(r, Q_cold, label="T=1e3 K", color='C0')

    plt.axvline(28.5, color="gray", ls="--")
    plt.axhline(1, color='gray', ls='--')

    plt.ylim(0.1,1e3)
    plt.xlim(0,30)

    plt.ylabel("Q")
    plt.xlabel("r [kpc]")

    plt.legend()
    plt.savefig("toomre.png")