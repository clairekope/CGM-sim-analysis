#######################################################
# Save disk masses (gas, stars, total) for all datasets
#######################################################

import glob
import yt
from yt.units import Myr
from yt.utilities.exceptions import YTFieldNotFound
import numpy as np
import calc_sfr
yt.enable_parallelism()

# Calculate SFR
latest = sorted(glob.glob("DD????/DD????"))[-1]

dt = 50 * Myr
ds_num = int(latest[-4:])
final_time = ds_num * dt# + dt/2
nbins = ds_num + 1 #2
time = np.linspace(0*Myr, final_time, nbins)

# Calculate masses over time
datasets = yt.load("DD????/DD????")
storage = {}
for my_storage, ds in datasets.piter(dynamic=False, storage=storage):

    ds.index # load field index

    if all([field in ds.field_list for field in calc_sfr.required_fields]):

        ds.add_field(name=('io','particle_initial_mass'),
                     units='g',
                     sampling_type='particle',
                     function=calc_sfr.field_initial_mass)

        ad = ds.all_data()
        formed_mass = ad.quantities.total_quantity(('io','particle_initial_mass'))
        star_mass = ad.quantities.total_quantity(('io','particle_mass'))

    # 1.3 kpc is 4 scale heights. Need more than that in radial.
    dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (20,'kpc'), (1.3, 'kpc'))
    disk_mass = dsk.quantities.total_quantity(('gas','cell_mass'))

    data = {}
    data['disk_mass'] = disk_mass.to('Msun')
    data['formed_mass'] = formed_mass.to('Msun')
    data['star_mass'] = star_mass.to('Msun')
    data['total_mass'] = (disk_mass + star_mass).to('Msun')
    
    my_storage.result = data
    my_storage.result_id = int(ds.basename[-4:])

# Compile and save data
if yt.is_root():
    data_arr = np.zeros((len(storage), 5))
    for ds_name, dic in storage.items():
        data_arr[ds_name, 0] = time[ds_name]
        data_arr[ds_name, 1] = dic['disk_mass']
        data_arr[ds_name, 2] = dic['formed_mass']
        data_arr[ds_name, 3] = dic['total_mass']
        data_arr[ds_name, 4] = dic['star_mass']

    np.savetxt('masses_over_time.txt', data_arr, header="Time_Myr DiskGas_Msun StellarMass_Msun Sum_Msun FormedMass_Msun")
