#!/usr/bin/env python
# coding: utf-8
import h5py
import napari
import numpy as np
import unyt as u
import scipy.stats as s
from scipy.interpolate import interp1d
from tqdm import tqdm

class Perturbations():

    def __init__(self, rng_seed=None, nsamples=500, 
                 grid_min=-800, grid_max=800, grid_size=512,
                 entropy_filename="../initial_entropy_profile.txt",
                 pressure_filename="../initial_pressure_profile.txt",
                 perturbation_filename="entropy_perturbations.h5"):
        
        self.rng_seed = rng_seed
        self.nsamples = nsamples
        self.grid_min = grid_min
        self.grid_max = grid_max
        self.grid_size = grid_size
        self.entropy_filename = entropy_filename
        self.pressure_filename = pressure_filename
        self.perturbation_filename = perturbation_filename

        self.rng = np.random.default_rng(self.rng_seed)

        self.entropy_func = None
        self.pressure_func = None

        self.x = np.linspace(self.grid_min, self.grid_max, self.grid_size)
        self.xx, self.yy, self.zz = np.meshgrid(self.x, self.x, self.x)
        self.rr = np.sqrt(self.xx**2 + self.yy**2 + self.zz**2)

        self.modified_entropy = None

    def load_entropy(self):

        if self.entropy_func is None:
            func_r, func_K = np.genfromtxt(self.entropy_filename, unpack=True)
            self.entropy_func = interp1d(func_r, func_K, fill_value="extrapolate")
    
    def load_pressure(self):

        if self.pressure_func is None:
            func_r, func_P = np.genfromtxt(self.pressure_filename, unpack=True)
            self.pressure_func = interp1d(func_r, func_P, fill_value="extrapolate")

    def generate_perturbations(self,
                               dlnK_sigma=0.6, dlnK_loc=-0.2, dlnK_scale=0.3,
                               tcross_mean=0.3, tcross_sigma=0.01):

        # Pull in the entropy and pressure profiles K(r) and P(r)
        self.load_entropy()
        self.load_pressure()


        # Randomly pick an (x,y,z) coordinate and a max dlnK for the perturbation
        # 
        # The use of a lognormal distribution for dlnK as well as the choice of sigma_lnK are based on Voit 2021:
        # 
        # >"Around each median profile, fluctuations in density,
        # temperature, and entropy exhibit approximately log-normal
        # distributions with long tails toward lower temperature and
        # entropy and greater gas density"
        # 

        pick_x = self.rng.choice(self.x, self.nsamples)
        pick_y = self.rng.choice(self.x, self.nsamples)
        pick_z = self.rng.choice(self.x, self.nsamples)

        pick_r = np.sqrt(pick_x**2 + pick_y**2 + pick_z**2)

        lognorm_params = {"s":dlnK_sigma, 
                        "loc":dlnK_loc,
                        "scale":dlnK_scale, 
                        "size":self.nsamples,
                        "random_state":self.rng}

        pick_dlnK = -1*s.lognorm.rvs(**lognorm_params)

        # Cut off at dlnK = -1
        while not np.all(pick_dlnK > -1): 
            pick_dlnK = np.where(pick_dlnK > -1, pick_dlnK, -1*s.lognorm.rvs(**lognorm_params))


        # Randomly pick a size for the perturbation using normally distributed sound crossing times.
        # The value we pick from that Gaussian distribution will be 2sigma of a Gaussian blur (roughly analogous to a diameter, but OK if underresolved)

        # assume n_e = n_i s.t. n = 2n_e
        gamma = 5/3
        mu = 0.6
        n = np.power(self.pressure_func(pick_r) / ( 2*self.entropy_func(pick_r) ), 1/gamma) / u.cm**3
        c_s = np.sqrt(gamma * self.pressure_func(pick_r)*u.erg/u.cm**3 / (2*mu*u.mp*n))

        pick_tcross = self.rng.normal(tcross_mean, tcross_sigma, self.nsamples) * u.Myr

        assert not (pick_tcross < 0).any()

        diam = (c_s * pick_tcross).to("kpc")

        # Save x,y,z dlnK, and simga (diam) to HDF5, because it's easier to get into enzo
        with h5py.File(self.perturbation_filename, "w") as f:
            f.create_dataset("x", data=pick_x)
            f.create_dataset("y", data=pick_y)
            f.create_dataset("z", data=pick_z)
            f.create_dataset("dlnK", data=pick_dlnK)
            f.create_dataset("sigma", data=diam)

    def apply_perturbations(self, overwrite=False):

        if self.modified_entropy is not None and overwrite is False:
            print("Modified entropy profile already saved to this object. Pass overwrite=True to recalculate.")
            return

        self.load_entropy()

        with h5py.File(self.perturbation_filename,"r") as f:
            x = np.array(f["x"])
            y = np.array(f["y"])
            z = np.array(f["z"])
            dlnK = np.array(f["dlnK"])
            sigma = np.array(f["sigma"])
        
        self.modified_entropy = self.entropy_func(self.rr)

        for i in tqdm(range(self.nsamples)):
            profile = dlnK[i] * np.exp(-(self.xx-x[i])**2/(2*(sigma[i]/2)**2) \
                                    -(self.yy-y[i])**2/(2*(sigma[i]/2)**2) \
                                    -(self.zz-z[i])**2/(2*(sigma[i]/2)**2) )
            self.modified_entropy = self.modified_entropy * profile + self.modified_entropy

    def save_perturbed_entropy(self, filename="entropy_perturbed.npy"):

        if self.modified_entropy is None:
            self.apply_perturbations()

        np.save(filename, self.modified_entropy, allow_pickle=False)

    def load_perturbed_entropy(self, filename="entropy_perturbed.npy", overwrite=False):

        if self.modified_entropy is not None and overwrite is False:
            print("Modified entropy profile already saved to this object. Pass overwrite=True to continue loading.")
            return

        temp = np.load(filename, allow_pickle=False)

        if temp.shape != (self.grid_size, self.grid_size, self.grid_size):
            raise RuntimeError(f"Loaded array with shape {temp.shape} does not match the cubic grid_size {self.grid_size} specified for this object")
        else:
            self.modified_entropy = temp

    def visualize(self, cmap="viridis"):
        
        if self.modified_entropy is None:
            print("This object has no modified entropy profile stored. Either run apply_perturbations() or load_perturbed_entropy().")
            return
        
        viewer = napari.view_image(self.modified_entropy, colormap=cmap)
        napari.run()

if __name__ == "__main__":
    
    p = Perturbations()
    p.generate_perturbations()
    p.save_perturbed_entropy()
    p.visualize()
