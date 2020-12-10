import glob
import yt
from yt.units import Myr
import numpy as np
from calc_sfr import calc_sfr
yt.enable_parallelism()

# Calculate SFR
latest = sorted(glob.glob("DD????/DD????"))[-1]

#ds = yt.load(latest)
#dt = ds.quan(50, 'Myr')
#ds_num = int(ds.basename[-4:])
#start_time = -dt/2

dt = 50 * Myr
ds_num = int(latest[-4:])
final_time = ds_num * dt# + dt/2
nbins = ds_num + 1 #2
time = np.linspace(0*Myr, final_time, nbins)

# contruct time bins so that bin centers match output times
#bins = np.linspace(start_time, final_time, nbins)
#time, sfr = calc_sfr(ds.all_data(), bins.to('yr'))
#time = ds.arr(time, 'yr').to('Myr')

# Calculate masses over time
datasets = yt.load("DD????/DD????")
storage = {}
for my_storage, ds in datasets.piter(storage=storage):

    # Field added by calc_sfr.py on import
    star_mass = ds.all_data().quantities.total_quantity('particle_initial_mass')
    
    # 4 scale heights
    dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (14,'kpc'), (1.3, 'kpc'))
    disk_mass = dsk.quantities.total_quantity('cell_mass')

    data = {}
    # data['time'] = ds.current_time.to('Gyr')
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
        data_arr[ds_name, 1] = data['disk_mass']
        data_arr[ds_name, 2] = data['star_mass']
        data_arr[ds_name, 3] = data['total_mass']

    np.savetxt('masses_over_time.txt', data_arr, header="Time_Myr DiskGas_Msun StellarMass_Msun Sum_Msun")
