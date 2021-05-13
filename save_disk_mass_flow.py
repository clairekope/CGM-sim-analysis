import yt
import numpy as np
yt.enable_parallelism()

datasets = yt.load("DD????/DD????")
storage = {}
for my_storage, ds in datasets.piter(dynamic=False, storage=storage):

    dsk = ds.disk([0.5,0.5,0.5], [0,0,1], (20,'kpc'), (1.3, 'kpc'))
    
    # 10 kpc wider
    buffer1 = ds.disk([0.5,0.5,0.5], [0,0,1], (30,'kpc'), (1.3, 'kpc'))
    shl = buffer1 - dsk
    
    # 5 kpc taller in each direction
    buffer2 = ds.disk([0.5,0.5,0.5], [0,0,1], (20,'kpc'), (6.3, 'kpc'))
    hat = buffer2 - dsk
    
    cold_dsk = dsk.cut_region(["obj['temperature'] < 1e4"])
    warm_dsk = dsk.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e6)"])
    hot_dsk = dsk.cut_region(["obj['temperature'] > 1e6"])
    
    cold_shl = shl.cut_region(["obj['temperature'] < 1e4"])
    warm_shl = shl.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e6)"])
    hot_shl = shl.cut_region(["obj['temperature'] > 1e6"])
    
    cold_hat = hat.cut_region(["obj['temperature'] < 1e4"])
    warm_hat = hat.cut_region(["(obj['temperature'] > 1e4) & (obj['temperature'] < 1e6)"])
    hot_hat = hat.cut_region(["obj['temperature'] > 1e6"])
    
    
    my_storage.result_id = int(ds.basename[-4:])
    my_storage.result = [cold_dsk['cell_mass'].sum(),
                         cold_shl['cell_mass'].sum(),
                         cold_hat['cell_mass'].sum(),
                         warm_dsk['cell_mass'].sum(),
                         warm_shl['cell_mass'].sum(),
                         warm_hat['cell_mass'].sum(),
                         hot_dsk['cell_mass'].sum(),
                         hot_shl['cell_mass'].sum(),
                         hot_hat['cell_mass'].sum()]
    
if yt.is_root():
    data_arr = np.zeros((len(storage), 9))
    for ds_name, lst in storage.items():
        print(lst)
        for i in range(9):
            data_arr[ds_name, i] = lst[i].to("Msun").v
            
    header = "ColdDiskMsun ColdShellMsun ColdHatMsun " + \
             "WarmDiskMsun WarmShellMsun WarmHatMsun " + \
             "HotDiskMsun HotShellMsun HotHatMsun"
    np.savetxt("disk_mass_flow.txt", data_arr, header=header)