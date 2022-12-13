#############################################
# Compile distributions of mass below a given
# cooling time for each dataset into a table
#############################################

import yt
yt.enable_parallelism()

import yt.units as u
import numpy as np

datasets = yt.load("DD????/DD????")
storage = {}

cgm_only = True

for my_storage, ds in datasets.piter(dynamic=False, storage=storage):
    sph = ds.sphere('c',(206,'kpc'))
       
    if cgm_only:
        dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (20,'kpc'), (1.3, 'kpc'))
        cgm = sph - dsk
        cgm.set_field_parameter('center', ds.arr([0.5,0.5,0.5], 'code_length'))
        dataset = cgm
    else:
        dataset = sph

    # Cummulative profile
    prof = yt.create_profile(dataset, ("gas","cooling_time"), ("gas","cell_mass"),
                             accumulation=False, fractional=False, weight_field=None,
                             units={'cooling_time':'Gyr','cell_mass':'Msun'},
                             extrema={"cooling_time":(1e-6, 1e6)})

    # Std Dev of tcool as function of tcool
    tc_binner = np.digitize(dataset['gas','cooling_time'], prof.x_bins)
    tcool_stddev = np.zeros(prof.x.size)
    for i in range(1, prof.x_bins.size):
        this_bin = tc_binner==i
        tcool_stddev[i-1] = np.std(dataset['gas','cooling_time'][this_bin].to('Gyr'))

    my_storage.result_id = int(ds.basename[-4:])
    my_storage.result = (prof['cell_mass'], tcool_stddev)
    
if yt.is_root():
    data_arr = np.zeros((prof.x.size, len(datasets)+1))
    stddev_arr = np.zeros_like(data_arr)

    data_arr[:,0] = prof.x
    stddev_arr[:,0] = prof.x
    
    for i in range(len(datasets)):
        data_arr[:,i+1] = storage[0][i]
        stddev_arr[:,i+1] = storage[1][i]
        
    np.savetxt("tcool_mass_dist_CGM{}.txt".format("" if cgm_only else "-disk"),
               data_arr,
               header="tcool [Gyr]  DD0000 mass [Msun]  DD0001 mass ... etc")

    np.savetxt("tcool_stddev_CGM{}.txt".format("" if cgm_only else "-disk"),
               stddev_arr,
               header="tcool [Gyr] DD0000 stddev(tcool) [Gyr] DD0001 stddev(tcool) ... etc")