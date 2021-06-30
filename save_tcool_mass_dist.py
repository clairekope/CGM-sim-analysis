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

for my_storage, ds in datasets.piter(dynamic=False, storage=storage):
    sph = ds.sphere('c',(206,'kpc'))
    prof = yt.create_profile(sph, "cooling_time", "cell_mass", 
                             accumulation=True,
                             units={'cooling_time':'Gyr','cell_mass':'Msun'},
                             extrema={"cooling_time":(1e-6, 1e6)})
    my_storage.result_id = int(ds.basename[-4:])
    my_storage.result = prof['cell_mass']
    
if yt.is_root():
    data_arr = np.zeros((prof.x.size, len(datasets)+1))
    data_arr[:,0] = prof.x
    for i in range(len(datasets)):
        data_arr[:,i+1] = storage[i]
        
    np.savetxt("tcool_mass_dist.txt", data_arr,
               header="tcool [Gyr]  DD0000 mass [Msun]  DD0001 mass ... etc")