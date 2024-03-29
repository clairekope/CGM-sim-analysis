##########################################################
# For determining where the disk gas & stars are,
# save narrow FOV density projections & the star positions
##########################################################

import yt
import numpy as np

for dataset in ["DD0006/DD0006",
                "DD0012/DD0012",
                "DD0018/DD0018",
                "DD0046/DD0046",
                "DD0052/DD0052",
                "DD0058/DD0058",
                "DD0066/DD0066",
                "DD0072/DD0072",
                "DD0078/DD0078"]:
    ds = yt.load(dataset)

    width_x = (80, "kpc")
    width_z = (80, "kpc")

    p_x = yt.ProjectionPlot(ds, 'x', 'density', width=width_x)
    p_z = yt.ProjectionPlot(ds, 'z', 'density', width=width_z)

    frb_x = p_x.data_source.to_frb(width_x, 512)
    frb_z = p_z.data_source.to_frb(width_z, 512)

    d_x = np.array(frb_x["density"])
    d_z = np.array(frb_z["density"])

    np.save(dataset[-6:]+"_x_density.npy", d_x)
    np.save(dataset[-6:]+"_z_density.npy", d_z)

    center = ds.quan(0.5, "code_length")
    star_coords = np.column_stack(((ds.r['particle_position_x']-center).to('kpc'),
                                   (ds.r['particle_position_y']-center).to('kpc'),
                                   (ds.r['particle_position_z']-center).to('kpc'),
                                   ds.r['particle_mass'].to('Msun')))
    np.savetxt(dataset[-6:]+"_star_coords.txt", star_coords,
               header="x y z")
