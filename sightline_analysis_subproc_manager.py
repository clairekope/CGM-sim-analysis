#################################################
# Manage processes running sightline_analysis.py,
# as it likes to read datasets wrong & segfault
#################################################

import yt
yt.enable_parallelism(suppress_logging=True)
import glob
import subprocess
import os
from mpi4py import MPI
from sightline_analysis import process_dataset

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def launch_subprocess(filename):
    proc = subprocess.Popen([#'ibrun',
                             #'-n',
                             #'1',
                             'python',
                             '/home/kopec/CGM-sim-analysis/sightline_analysis.py',
                             filename])
    return proc


filenames = glob.glob("DD??[0-8]?/DD????")
for filename in yt.parallel_objects(filenames):

    #print("Parent CPU:", os.sched_getaffinity(0))
    
    rays_bad = True

    while rays_bad:
        proc = launch_subprocess(filename)
        print("Rank", rank, "launched", filename, "on pid", proc.pid, flush=True)

        proc.wait() # wait for subprocces to finish

        print("Rank", rank, "pip", proc.pid, "returncode:", proc.returncode, flush=True)
        rays_bad = proc.returncode # anything but 0 is failure
        
    print("Rank", rank, "finished file", filename, flush=True)
