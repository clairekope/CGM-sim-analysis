from blk import pipeline, cache
from perturbation_generator import Perturbations

cachedir = "blk_test/"

cache.set_dir("blk_test/")

p = Perturbations(rng_seed=333, nsamples=500, 
                  perturbation_filename=cachedir+"entropy_perturbations.h5")

gen = pipeline.Stage(p.generate_perturbations(),
                     arguments=[],
                     tag = "generation",
                     action = "manual")

app = pipeline.Stage(p.save_perturbed_entropy(),
                     arguments=[cachedir+"entropy_perturbed.npy"],
                     tag = "application",
                     action = "manual",
                     depends_on = [gen])

pipeline.execute(app, parallelism="none")

p.load_perturbed_entropy(cachedir+"entropy_perturbed.npy")
p.visualize()