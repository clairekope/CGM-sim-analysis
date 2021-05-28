import yt
from yt import derived_field
import numpy as np
yt.enable_parallelism()
import pdb

@derived_field(name=("gas","momentum_cylindrical_z"), units="g*cm/s")
def _vertical_flow(field, data):
    return data['cell_mass'] * data['velocity_cylindrical_z']

@derived_field(name=("gas","momentum_cylindrical_radius"), units="g*cm/s")
def _radial_flow(field, data):
    return data['cell_mass'] * data['velocity_cylindrical_radius']

datasets = yt.load("DD????/DD????", unit_system='galactic')
storage = {}

field_list = ['density','momentum_cylindrical_z','momentum_cylindrical_radius']

for my_storage, ds in datasets.piter(dynamic=False, storage=storage):

    dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (20,'kpc'), (1.3, 'kpc'), fields=field_list)
    
    # 10 kpc wider
    buffer1 = ds.disk([0.5,0.5,0.5], [0,0,1], (30,'kpc'), (1.3, 'kpc'), fields=field_list)
    rng = buffer1 - dsk
    
    # 5 kpc taller in each direction
    buffer2 = ds.disk([0.5,0.5,0.5], [0,0,1], (20,'kpc'), (6.3, 'kpc'), fields=field_list)
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
    my_storage.result = [cold_dsk.quantities.total_quantity('cell_mass'),
                         cold_rng.quantities.total_quantity('cell_mass'),
                         cold_slb.quantities.total_quantity('cell_mass'),
                         warm_dsk.quantities.total_quantity('cell_mass'),
                         warm_rng.quantities.total_quantity('cell_mass'),
                         warm_slb.quantities.total_quantity('cell_mass'),
                         hot_dsk.quantities.total_quantity('cell_mass'),
                         hot_rng.quantities.total_quantity('cell_mass'),
                         hot_slb.quantities.total_quantity('cell_mass'),
                         cold_dsk.quantities.total_quantity('momentum_cylindrical_z'),
                         cold_rng.quantities.total_quantity('momentum_cylindrical_z'),
                         cold_slb.quantities.total_quantity('momentum_cylindrical_z'),
                         warm_dsk.quantities.total_quantity('momentum_cylindrical_z'),
                         warm_rng.quantities.total_quantity('momentum_cylindrical_z'),
                         warm_slb.quantities.total_quantity('momentum_cylindrical_z'),
                         hot_dsk.quantities.total_quantity('momentum_cylindrical_z'),
                         hot_rng.quantities.total_quantity('momentum_cylindrical_z'),
                         hot_slb.quantities.total_quantity('momentum_cylindrical_z'),
                         cold_dsk.quantities.total_quantity('momentum_cylindrical_radius'),
                         cold_rng.quantities.total_quantity('momentum_cylindrical_radius'),
                         cold_slb.quantities.total_quantity('momentum_cylindrical_radius'),
                         warm_dsk.quantities.total_quantity('momentum_cylindrical_radius'),
                         warm_rng.quantities.total_quantity('momentum_cylindrical_radius'),
                         warm_slb.quantities.total_quantity('momentum_cylindrical_radius'),
                         hot_dsk.quantities.total_quantity('momentum_cylindrical_radius'),
                         hot_rng.quantities.total_quantity('momentum_cylindrical_radius'),
                         hot_slb.quantities.total_quantity('momentum_cylindrical_radius'),]

    pdb.set_trace()
    
if yt.is_root():
    data_arr = np.zeros((len(storage), 27))
    for ds_name, lst in storage.items():
        for i in range(data_arr.shape[1]):
            data_arr[ds_name, i] = lst[i].v
            
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
