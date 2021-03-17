import glob
import yt
from yt.units import Myr
from yt.utilities.exceptions import YTFieldNotFound
import numpy as np
from calc_sfr import add_initial_mass_field
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
for my_storage, ds in datasets.piter(dynamic=True, storage=storage):

    try:
        add_initial_mass_field(ds)
        star_mass = ds.all_data().quantities.total_quantity('particle_initial_mass')
    except YTFieldNotFound:
        star_mass = ds.quan(0, 'Msun')

    # 4 scale heights
    dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (14,'kpc'), (1.3, 'kpc'))
    disk_mass = dsk.quantities.total_quantity('cell_mass')

    data = {}
    data['disk_mass'] = disk_mass.to('Msun')
    data['star_mass'] = star_mass.to('Msun')
    data['total_mass'] = (disk_mass + star_mass).to('Msun')
    
    my_storage.result = data
    my_storage.result_id = int(ds.basename[-4:])

# Compile and save data
if yt.is_root():
    data_arr = np.zeros((len(storage), 4))
    for ds_name, dic in storage.items():
        data_arr[ds_name, 0] = time[ds_name]
        data_arr[ds_name, 1] = dic['disk_mass']
        data_arr[ds_name, 2] = dic['star_mass']
        data_arr[ds_name, 3] = dic['total_mass']

    np.savetxt('masses_over_time.txt', data_arr, header="Time_Myr DiskGas_Msun StellarMass_Msun Sum_Msun")
