import yt
from yt import derived_field
import numpy as np
yt.enable_parallelism()

@derived_field(name=("gas","momentum_cylindrical_z"), units="g*cm/s")
def _vertical_flow(field, data):
    return data['cell_mass'] * data['velocity_cylindrical_z']

@derived_field(name=("gas","momentum_cylindrical_radius"), units="g*cm/s")
def _radial_flow(field, data):
    return data['cell_mass'] * data['velocity_cylindrical_radius']

datasets = yt.load("DD????/DD????")
storage = {}
for my_storage, ds in datasets.piter(dynamic=False, storage=storage):

    dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (20,'kpc'), (1.3, 'kpc'))
    
    # 10 kpc wider
    buffer1 = ds.disk([0.5,0.5,0.5], [0,0,1], (30,'kpc'), (1.3, 'kpc'))
    rng = buffer1 - dsk
    
    # 5 kpc taller in each direction
    buffer2 = ds.disk([0.5,0.5,0.5], [0,0,1], (20,'kpc'), (6.3, 'kpc'))
    slb = buffer2 - dsk
    
    cold_dsk = dsk.cut_region(["obj['temperature'] < 1e4"])
    warm_dsk = dsk.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e6)"])
    hot_dsk = dsk.cut_region(["obj['temperature'] > 1e6"])
    
    cold_rng = rng.cut_region(["obj['temperature'] < 1e4"])
    warm_rng = rng.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e6)"])
    hot_rng = rng.cut_region(["obj['temperature'] > 1e6"])
    
    cold_slb = slb.cut_region(["obj['temperature'] < 1e4"])
    warm_slb = slb.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e6)"])
    hot_slb = slb.cut_region(["obj['temperature'] > 1e6"])
    
#     cold_rad = buffer1.cut_region(["obj['temperature'] < 1e4"])
#     warm_rad = buffer1.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e6)"])
#     hot_rad = buffer1.cut_region(["obj['temperature'] > 1e6"])
    
#     cold_vrt = buffer2.cut_region(["obj['temperature'] < 1e4"])
#     warm_vrt = buffer2.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e6)"])
#     hot_vrt = buffer2.cut_region(["obj['temperature'] > 1e6"])
    
    my_storage.result_id = int(ds.basename[-4:])
    my_storage.result = [cold_dsk['cell_mass'].sum(),
                         cold_rng['cell_mass'].sum(),
                         cold_slb['cell_mass'].sum(),
                         warm_dsk['cell_mass'].sum(),
                         warm_rng['cell_mass'].sum(),
                         warm_slb['cell_mass'].sum(),
                         hot_dsk['cell_mass'].sum(),
                         hot_rng['cell_mass'].sum(),
                         hot_slb['cell_mass'].sum(),
                         cold_dsk[('gas','momentum_cylindrical_z')].sum(),
                         cold_rng[('gas','momentum_cylindrical_z')].sum(),
                         cold_slb[('gas','momentum_cylindrical_z')].sum(),
                         warm_dsk[('gas','momentum_cylindrical_z')].sum(),
                         warm_rng[('gas','momentum_cylindrical_z')].sum(),
                         warm_slb[('gas','momentum_cylindrical_z')].sum(),
                         hot_dsk[('gas','momentum_cylindrical_z')].sum(),
                         hot_rng[('gas','momentum_cylindrical_z')].sum(),
                         hot_slb[('gas','momentum_cylindrical_z')].sum(),
                         cold_dsk[('gas','momentum_cylindrical_radius')].sum(),
                         cold_rng[('gas','momentum_cylindrical_radius')].sum(),
                         cold_slb[('gas','momentum_cylindrical_radius')].sum(),
                         warm_dsk[('gas','momentum_cylindrical_radius')].sum(),
                         warm_rng[('gas','momentum_cylindrical_radius')].sum(),
                         warm_slb[('gas','momentum_cylindrical_radius')].sum(),
                         hot_dsk[('gas','momentum_cylindrical_radius')].sum(),
                         hot_rng[('gas','momentum_cylindrical_radius')].sum(),
                         hot_slb[('gas','momentum_cylindrical_radius')].sum(),]
    
if yt.is_root():
    data_arr = np.zeros((len(storage), 27))
    for ds_name, lst in storage.items():
        print(lst)
        for i in range(9):
            data_arr[ds_name, i] = lst[i].to("Msun").v
            
    header = "Disk_Cold_Mgas Ring_Cold_Mgas Slab_Cold_Mgas " + \
             "Disk_Warm_Mgas Ring_Warm_Mgas Slab_Warm_Mgas " + \
             "Disk_Hot_Mgas  Ring_Hot_Mgas  Slab_Hot_Mgas " + \
             "Disk_Cold_pz Ring_Cold_pz Slab_Cold_pz " + \
             "Disk_Warm_pz Ring_Warm_pz Slab_Warm_pz " + \
             "Disk_Hot_pz  Ring_Hot_pz  Slab_Hot_pz " + \
             "Disk_Cold_pr Ring_Cold_pr Slab_Cold_pr " + \
             "Disk_Warm_pr Ring_Warm_pr Slab_Warm_pr " + \
             "Disk_Hot_pr  Ring_Hot_pr  Slab_Hot_pr"
    np.savetxt("disk_mass_flow.txt", data_arr, header=header)