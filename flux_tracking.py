"""
Filename: flux_tracking.py
Author: Cassi
Date created: 9-27-19
Date last modified: 4-30-21
This file takes command line arguments and computes fluxes of things through surfaces.

Dependencies:
utils/consistency.py
utils/get_refine_box.py
utils/get_halo_center.py
utils/get_proper_box_size.py
utils/get_run_loc_etc.py
utils/yt_fields.py
utils/foggie_load.py
utils/analysis_utils.py
"""

import numpy as np
import yt
import unyt as u
from yt import YTArray, derived_field
import glob
from astropy.table import Table
import datetime
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from flux_utils import *
from calc_enclosed_mass import NFW_mass_enclosed

@derived_field(name="radial_kinetic_energy", sampling_type="cell", units="erg")
def _radial_kinetic_energy(field, data):
    return 0.5 * data['cell_mass'] * data['radial_velocity']**2.

@derived_field(name="tangential_kinetic_energy", sampling_type="cell", units="erg")
def _tangential_kinetic_energy(field, data):
    return 0.5 * data['cell_mass'] * data['tangential_velocity']**2.

def set_table_units(table):
    '''Sets the units for the table. Note this needs to be updated whenever something is added to
    the table. Returns the table.'''

    for key in table.keys():
        if ('radius' in key) or ('height' in key) or ('edge' in key):
            table[key].unit = 'kpc'
        elif ('mass' in key) or ('metal' in key):
            table[key].unit = 'Msun/yr'
        elif ('energy' in key):
            table[key].unit = 'erg/yr'
        elif ('entropy' in key):
            table[key].unit = 'cm**2*keV/yr'
        elif ('O' in key):
            table[key].unit = 'Msun/yr'
        elif ('momentum' in key):
            table[key].unit = 'g*cm**2/s/yr'
    return table

def make_table(flux_types, surface_type, temp_cut=False, edge=False):
    '''Makes the giant table that will be saved to file.
    If 'edge' is True, this table is for fluxes into/out of edges of shape.'''

    if (surface_type[0]=='sphere'):
        if (edge):
            names_list = ['inner_radius', 'outer_radius']
            types_list = ['f8', 'f8']
        else:
            names_list = ['radius']
            types_list = ['f8']
    if (surface_type[0]=='cylinder'):
        if (surface_type[1]=='radius'):
            if (edge):
                names_list = ['inner_radius', 'outer_radius']
                types_list = ['f8', 'f8']
            else:
                names_list = ['radius']
                types_list = ['f8']
        if (surface_type[1]=='height'):
            if (edge):
                names_list = ['bottom_edge', 'top_edge']
                types_list = ['f8', 'f8']
            else:
                names_list = ['height']
                types_list = ['f8']

    dir_name = ['net_', '_in', '_out']
    if (temp_cut): temp_name = ['', 'cold_', 'cool_', 'warm_', 'hot_']
    else: temp_name = ['']

    for i in range(len(flux_types)):
        for k in range(len(temp_name)):
            for j in range(len(dir_name)):
                if (flux_types[i]=='cooling_energy_flux') and (not edge):
                    if (j==0):
                        name = 'net_' + temp_name[k] + 'cooling_energy_flux'
                        names_list += [name]
                        types_list += ['f8']
                elif (flux_types[i]!='cooling_energy_flux'):
                    if (j==0): name = dir_name[j]
                    else: name = ''
                    name += temp_name[k]
                    name += flux_types[i]
                    if (j>0): name += dir_name[j]
                    names_list += [name]
                    types_list += ['f8']

    table = Table(names=names_list, dtype=types_list)

    return table

def calc_fluxes(ds, dt, tablename, save_suffix, 
                surface_args, flux_types, Menc_profile, disk=False,
                temp_cut=False, units_rvir=False, Rvir=100.):
    '''This function calculates the fluxes specified by 'flux_types' into and out of the surfaces specified by 'surface_args'. 
    It uses the dataset stored in 'ds' and the time step between outputs is 'dt'. 
    Fluxes are stored in 'tablename' with 'save_suffix' appended. If 'disk' is True,
    then at least once surface shape requires disk-relative fields.

    This function calculates the flux as the sum of all cells whose velocity and distance from the
    surface of interest indicate that the gas contained in that cell will be displaced across the
    surface of interest by the next timestep. That is, the properties of a cell contribute to the
    flux if it is no further from the surface of interest than v*dt where v is the cell's velocity
    normal to the surface and dt is the time between snapshots, which is dt = 5.38e6 yrs for the DD
    outputs.

    This function differs from calc_fluxes below in that it returns just one flux across each surface
    specified, rather than fluxes at a number of steps within the total shape bounded by the surface.
    It only tracks things entering or leaving the closed surface and nothing else.'''

    # Set up table of everything we want
    fluxes = []
    flux_filename = ''
    if ('mass' in flux_types):
        fluxes.append('mass_flux')
        fluxes.append('metal_flux')
        flux_filename += '_mass'
    if ('energy' in flux_types):
        fluxes.append('thermal_energy_flux')
        fluxes.append('kinetic_energy_flux')
        fluxes.append('radial_kinetic_energy_flux')
        fluxes.append('tangential_kinetic_energy_flux')
        fluxes.append('potential_energy_flux')
        fluxes.append('bernoulli_energy_flux')
        fluxes.append('cooling_energy_flux')
        flux_filename += '_energy'
    if ('entropy' in flux_types):
        fluxes.append('entropy_flux')
        flux_filename += '_entropy'

# Define list of ways to chunk up the shape over radius or height
    edges = False
    if (surface_args[0][0]=='cylinder'):
        table = make_table(fluxes, ['cylinder', surface_args[0][7]], temp_cut=temp_cut)
        table_edge = make_table(fluxes, ['cylinder', surface_args[0][7]], edge=True, temp_cut=temp_cut)
        edges = True
        bottom_edge = surface_args[0][1]
        top_edge = surface_args[0][2]
        cyl_radius = surface_args[0][6]
        num_steps = surface_args[0][3]
        
        if (surface_args[0][7]=='height'):
            if (units_rvir):
                dz = (top_edge-bottom_edge)/num_steps*Rvir
                chunks = ds.arr(np.arange(bottom_edge*Rvir,top_edge*Rvir+dz,dz), 'kpc')
            else:
                dz = (top_edge-bottom_edge)/num_steps
                chunks = ds.arr(np.arange(bottom_edge,top_edge+dz,dz), 'kpc')
        else:
            if (units_rvir):
                dr = cyl_radius/num_steps*Rvir
                chunks = ds.arr(np.arange(0.,cyl_radius*Rvir+dr,dr), 'kpc')
            else:
                dr = cyl_radius/num_steps
                chunks = ds.arr(np.arange(0.,cyl_radius+dr,dr), 'kpc')
    else:
        inner_radius = surface_args[0][1]
        outer_radius = surface_args[0][2]
        num_steps = surface_args[0][3]
        table = make_table(fluxes, ['sphere', 0], temp_cut=temp_cut)
        
        if (surface_args[0][0]=='frustum') or (surface_args[0][0]=='ellipse'):
            table_edge = make_table(fluxes, ['sphere', 0], edge=True, temp_cut=temp_cut)
            edges = True

        if (units_rvir):
            dr = (outer_radius-inner_radius)/num_steps*Rvir
            chunks = ds.arr(np.arange(inner_radius*Rvir,outer_radius*Rvir+dr,dr), 'kpc')
        else:
            dr = (outer_radius-inner_radius)/num_steps
            chunks = ds.arr(np.arange(inner_radius,outer_radius+dr,dr), 'kpc')

    # Load arrays of all fields we need
    print('Loading field arrays')
    sphere = ds.sphere('c', chunks[-1])
    # Filter CGM?

    radius = sphere['index','radius'].in_units('kpc').v
    x = sphere['gas','x'].in_units('kpc').v
    y = sphere['gas','y'].in_units('kpc').v
    z = sphere['gas','z'].in_units('kpc').v
    if (disk):
        x_disk = sphere['gas','x_disk'].in_units('kpc').v
        y_disk = sphere['gas','y_disk'].in_units('kpc').v
        z_disk = sphere['gas','z_disk'].in_units('kpc').v
        vx_disk = sphere['gas','vx_disk'].in_units('km/s').v
        vy_disk = sphere['gas','vy_disk'].in_units('km/s').v
        vz_disk = sphere['gas','vz_disk'].in_units('km/s').v
        new_x_disk = x_disk + vx_disk*dt*(100./cmtopc*stoyr)
        new_y_disk = y_disk + vy_disk*dt*(100./cmtopc*stoyr)
        new_z_disk = z_disk + vz_disk*dt*(100./cmtopc*stoyr)
    theta = sphere['index','spherical_theta'].v
    phi = sphere['index','spherical_phi'].v
    vx = sphere['gas','velocity_x'].in_units('km/s').v
    vy = sphere['gas','velocity_y'].in_units('km/s').v
    vz = sphere['gas','velocity_z'].in_units('km/s').v
    rad_vel = sphere['gas','radial_velocity'].in_units('km/s').v
    new_x = x + vx*dt*(100*u.pc/u.yr).to('cm/s').v
    new_y = y + vy*dt*(100*u.pc/u.yr).to('cm/s').v
    new_z = z + vz*dt*(100*u.pc/u.yr).to('cm/s').v
    new_radius = np.sqrt(new_x**2. + new_y**2. + new_z**2.)
    new_theta = np.arccos(new_z/new_radius)
    new_phi = np.arctan2(new_y, new_x)
    temperature = np.log10(sphere['gas','temperature'].in_units('K').v)
    fields = []
    if ('mass' in flux_types):
        mass = sphere['gas','cell_mass'].in_units('Msun').v
        metal_mass = sphere['gas','metal_mass'].in_units('Msun').v
        fields.append(mass)
        fields.append(metal_mass)
    if ('energy' in flux_types):
        kinetic_energy = (sphere['gas','kinetic_energy_density']*sphere['index','cell_volume']).in_units('erg').v
        radial_kinetic_energy = sphere['gas','radial_kinetic_energy'].in_units('erg').v
        if (disk): tangential_kinetic_energy = sphere['gas','tangential_kinetic_energy_disk'].in_units('erg').v
        else: tangential_kinetic_energy = sphere['gas','tangential_kinetic_energy'].in_units('erg').v
        thermal_energy = (sphere['gas','cell_mass']*sphere['gas','thermal_energy']).in_units('erg').v
        potential_energy = -(u.physical_constants.G * Menc_profile(radius*u.kpc) / (radius*u.kpc) * sphere['gas','cell_mass']).in_units('erg').v
        bernoulli_energy = kinetic_energy + 5./3.*thermal_energy + potential_energy
        cooling_energy = thermal_energy/sphere['gas','cooling_time'].in_units('yr').v
        fields.append(thermal_energy)
        fields.append(kinetic_energy)
        fields.append(radial_kinetic_energy)
        fields.append(tangential_kinetic_energy)
        fields.append(potential_energy)
        fields.append(bernoulli_energy)
        fields.append(cooling_energy)
    if ('entropy' in flux_types):
        entropy = sphere['gas','entropy'].in_units('keV*cm**2').v
        fields.append(entropy)
    
    # Cut to just the shapes specified
    if (disk):
        if (surface_args[0][0]=='cylinder'):
            bool_inshapes, radius = segment_region(x, y, z, theta, phi, radius, surface_args,
                                                   x_disk=x_disk, y_disk=y_disk, z_disk=z_disk,
                                                   Rvir=Rvir, units_rvir=units_rvir)
            bool_inshapes_new, new_radius = segment_region(new_x, new_y, new_z, new_theta, new_phi, new_radius, surface_args,
                                                   x_disk=new_x_disk, y_disk=new_y_disk, z_disk=new_z_disk,
                                                   Rvir=Rvir, units_rvir=units_rvir)
        else:
            bool_inshapes = segment_region(x, y, z, theta, phi, radius, surface_args,
                                           x_disk=x_disk, y_disk=y_disk, z_disk=z_disk, 
                                           Rvir=Rvir, units_rvir=units_rvir)
            bool_inshapes_new = segment_region(new_x, new_y, new_z, new_theta, new_phi, new_radius, surface_args,
                                           x_disk=new_x_disk, y_disk=new_y_disk, z_disk=new_z_disk, 
                                           Rvir=Rvir, units_rvir=units_rvir)
    else:
        if (surface_args[0][0]=='cylinder'):
            bool_inshapes, radius = segment_region(x, y, z, theta, phi, radius, surface_args,
                                                   Rvir=Rvir, units_rvir=units_rvir)
            bool_inshapes_new, new_radius = segment_region(new_x, new_y, new_z, new_theta, new_phi, new_radius, surface_args,
                                                   Rvir=Rvir, units_rvir=units_rvir)
        else:
            bool_inshapes = segment_region(x, y, z, theta, phi, radius, surface_args,
                                           Rvir=Rvir, units_rvir=units_rvir)
            bool_inshapes_new = segment_region(new_x, new_y, new_z, new_theta, new_phi, new_radius, surface_args,
                                               Rvir=Rvir, units_rvir=units_rvir)
    bool_inshapes_entire = (bool_inshapes) & (bool_inshapes_new)
    bool_toshapes = (~bool_inshapes) & (bool_inshapes_new)
    bool_fromshapes = (bool_inshapes) & (~bool_inshapes_new)

    # Cut to within shapes and entering/leaving shapes
    fields_shapes = []
    if (edges):
        fields_in_shapes = []
        fields_out_shapes = []
    radius_shapes = radius[bool_inshapes_entire]
    new_radius_shapes = new_radius[bool_inshapes_entire]
    temperature_shapes = temperature[bool_inshapes_entire]
    if (edges):
        radius_in_shapes = radius[bool_toshapes]
        new_radius_in_shapes = new_radius[bool_toshapes]
        temperature_in_shapes = temperature[bool_toshapes]
        radius_out_shapes = radius[bool_fromshapes]
        new_radius_out_shapes = new_radius[bool_fromshapes]
        temperature_out_shapes = temperature[bool_fromshapes]
    for i in range(len(fields)):
        field = fields[i]
        fields_shapes.append(field[bool_inshapes_entire])
        if (edges):
            fields_in_shapes.append(field[bool_toshapes])
            fields_out_shapes.append(field[bool_fromshapes])

    if (temp_cut): temps = [0.,4.,5.,6.,12.]
    else: temps = [0.]

    # Loop over chunks and compute fluxes to add to tables
    for r in range(len(chunks)-1):

        inner = chunks[r]
        outer = chunks[r+1]
        row = [inner]
        if (edges): row_edge = [inner, outer]

        temp_up = temperature_shapes[(radius_shapes < inner) & (new_radius_shapes > inner)]
        temp_down = temperature_shapes[(radius_shapes > inner) & (new_radius_shapes < inner)]
        temp_r = temperature_shapes[(radius_shapes > inner) & (radius_shapes < outer)]
        if (edges):
            temp_in = temperature_in_shapes[(new_radius_in_shapes > inner) & (new_radius_in_shapes < outer)]
            temp_out = temperature_out_shapes[(radius_out_shapes > inner) & (radius_out_shapes < outer)]

        for i in range(len(fields)):
            if (fluxes[i]=='cooling_energy_flux'):
                iter = [0]
                field_r = fields_shapes[i][(radius_shapes > inner) & (radius_shapes < outer)]
            else:
                iter = [0,1,2]
                field_up = fields_shapes[i][(radius_shapes < inner) & (new_radius_shapes > inner)]
                field_down = fields_shapes[i][(radius_shapes > inner) & (new_radius_shapes < inner)]
                if (edges):
                    field_in = fields_in_shapes[i][(new_radius_in_shapes > inner) & (new_radius_in_shapes < outer)]
                    field_out = fields_out_shapes[i][(radius_out_shapes > inner) & (radius_out_shapes < outer)]

            for k in range(len(temps)):
                if (k==0):
                    if (fluxes[i]=='cooling_energy_flux'):
                        field_r_t = field_r
                    else:
                        field_up_t = field_up
                        field_down_t = field_down
                        if (edges):
                            field_in_t = field_in
                            field_out_t = field_out
                else:
                    if (fluxes[i]=='cooling_energy_flux'):
                        field_r_t = field_r[(temp_r > temps[k-1]) & (temp_r < temps[k])]
                    else:
                        field_up_t = field_up[(temp_up > temps[k-1]) & (temp_up < temps[k])]
                        field_down_t = field_down[(temp_down > temps[k-1]) & (temp_down < temps[k])]
                        if (edges):
                            field_in_t = field_in[(temp_in > temps[k-1]) & (temp_in < temps[k])]
                            field_out_t = field_out[(temp_out > temps[k-1]) & (temp_out < temps[k])]
                for j in iter:
                    if (j==0):
                        if (fluxes[i]=='cooling_energy_flux'):
                            row.append(-np.sum(field_r_t))
                        else:
                            row.append(np.sum(field_up_t)/dt - np.sum(field_down_t)/dt)
                            if (edges):
                                row_edge.append(np.sum(field_in_t)/dt - np.sum(field_out_t)/dt)
                    if (j==1):
                        row.append(-np.sum(field_down_t)/dt)
                        if (edges):
                            row_edge.append(np.sum(field_in_t)/dt)
                    if (j==2):
                        row.append(np.sum(field_up_t)/dt)
                        if (edges):
                            row_edge.append(-np.sum(field_out_t)/dt)


        table.add_row(row)
        if (edges): table_edge.add_row(row_edge)

    table = set_table_units(table)
    if (edges): table_edge = set_table_units(table_edge)

    # Save to file
    table.write(tablename + flux_filename + save_suffix + '.hdf5', path='all_data', 
                serialize_meta=True, overwrite=True)
    if (edges): table_edge.write(tablename + '_edge' + flux_filename + save_suffix + '.hdf5', 
                                 path='all_data', serialize_meta=True, overwrite=True)

    return "Fluxes have been calculated for " + ds.basename

def load_and_calculate(tablename, save_suffix, surface_args, flux_types, 
                       temp_cut=False, units_rvir=False, Rvir=100.):
    '''This function loads a specified snapshot 'snap' located in the 'run_dir' within the
    'foggie_dir', the halo track 'track', the name of the halo_c_v file, the name of the snapshot,
    the name of the table to output, the mass enclosed table, the list of surface arguments, and
    the directory where the satellites file is saved, then
    does the calculation on the loaded snapshot.'''


    # Load the snapshot depending on if disk minor axis is needed
    disk = False
    for i in range(len(surface_args)):
        if (((surface_args[i][0]=='frustum') or (surface_args[i][0]=='cylinder')) and (surface_args[i][4]=='disk minor axis')):
            disk = True
    # load dataset
    if (disk):
        pass
        # define disk object

    dt = 5e7 # 50 Myr

    # Do the actual calculation
    message = calc_fluxes(ds, dt, tablename, save_suffix,
                                 surface_args, flux_types, NFW_mass_enclosed, disk=disk,
                                 temp_cut=temp_cut, units_rvir=units_rvir, Rvir=Rvir)

    print(message)
    print(str(datetime.datetime.now()))


if __name__ == "__main__":
    yt.enable_parallelism()

    surface = "['sphere',6,206,100]"
    flux_type = "energy,mass,entropy"
    save_suffix = ''
    units_rvir = False
    temp_cut = True

    surfaces = identify_shape(surface, units_rvir=units_rvir)

    # Build flux type list
    if (',' in flux_type):
        flux_types = flux_type.split(',')
    else:
        flux_types = [flux_type]
    for i in range(len(flux_types)):
        if (flux_types[i]!='mass') and (flux_types[i]!='energy') and (flux_types[i]!='entropy'):
            raise RuntimeError('The flux type   %s   has not been implemented. Ask Cassi to add it.' % (flux_types[i]))

    # Loop over outputs, for either single-processor or parallel processor computing
    datasets = yt.load('DD????/DD????')
    for ds in datasets.piter():
        # Make the output table name for this snapshot
        tablename = 'fluxes_'+ds.basename
        # Do the actual calculation
        load_and_calculate(tablename, save_suffix, surfaces, flux_types, temp_cut, units_rvir)

    print(str(datetime.datetime.now()))
    print("All snapshots finished!")
