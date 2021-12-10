import yt
import numpy as np
yt.enable_parallelism()

datasets = yt.load("DD????/DD????")
storage = {}
for my_storage, ds in datasets.piter(dynamic=False, storage=storage):

    streams = ds.cut_region(ds.all_data(), ["(obj['radial_velocity']).in_units('km/s') < 100"])
    cold = streams.include_below("temperature", 1e4)

    slc = yt.SlicePlot(ds, 'x', 'radial_velocity', width=(300,'kpc'), data_source=cold)
    slc.set_log('radial_velocity', False)
    slc.set_unit('radial_velocity', 'km/s')
    slc.set_zlim('radial_velocity',-60, 60)
    slc.set_cmap('radial_velocity', 'seismic')
    slc.save()

    var_x, mean_x = cold.quantities.weighted_variance('velocity_x', weight='cell_mass').in_units('km/s')
    var_y, mean_y = cold.quantities.weighted_variance('velocity_y', weight='cell_mass').in_units('km/s')
    var_z, mean_z = cold.quantities.weighted_variance('velocity_z', weight='cell_mass').in_units('km/s')

    vel_disp = var_x+var_y+var_z

    var_r, mean_r = cold.quantities.weighted_variance('radial_velocity', weight='cell_mass').in_units('km/s')

    data = {}
    data['vel_disp_3D'] = vel_disp
    data['vel_disp_radial'] = var_r
    data['mean_x'] = mean_x
    data['mean_y'] = mean_y
    data['mean_z'] = mean_z
    data['mean_r'] = mean_r
    
    my_storage.result = data
    my_storage.result_id = int(ds.basename[-4:])
    
if yt.is_root():
    data_arr = np.zeros((len(storage), 6))
    for ds_name, dic in storage.items():
        data_arr[ds_name, 0] = dic['vel_disp_3D']
        data_arr[ds_name, 0] = dic['vel_disp_r']
        data_arr[ds_name, 0] = dic['mean_x']
        data_arr[ds_name, 0] = dic['mean_y']
        data_arr[ds_name, 0] = dic['mean_y']
        data_arr[ds_name, 0] = dic['mean_y']
        
    np.savetxt("cold_gas_vel_disp.txt", data_arr, header="Disp3D DispR MeanX MeanY MeanZ MeanR")
    