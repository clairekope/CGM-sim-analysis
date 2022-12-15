"""
Filename: analysis_utils.py
This file contains functions that are used across all of Cassi's code, which includes:
-flux_tracking/flux_tracking.py
-radial_quantities/stats_in_shells.py
-radial_quantities/totals_in_shells.py
-plots/plot_fluxes.py
-plots/plot_1Dhistograms.py
-plots/plot_2Dhistograms.py
-plots/plot_outputs.py
-segmenting_regions/find_shape_for_region.py
-segmenting_regions/stack_FRBs.py
-paper_plots/mod_vir_temp/mod_vir_temp_paper_plots.py
Please let Cassi know if you make changes to this file!!!!!
"""

# Import everything as needed
from __future__ import print_function

import numpy as np
import yt
from yt.units import *
from yt import YTArray
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import ast

def identify_shape(shape_args, units_rvir=False):
    '''Returns an organized list of shape arguments from the input shape args.'''

    try:
        shape_args = ast.literal_eval(shape_args)
    except ValueError:
        sys.exit("Something's wrong with your shape arguments. Make sure to include both the outer " + \
        "quotes and the inner quotes around the shape type, like so:\n" + \
        '"[\'sphere\', 0.05, 2., 200.]"')
    if (type(shape_args[0])==str):
        shape_args = [shape_args]
    shapes = []
    for i in range(len(shape_args)):

        if (shape_args[i][0]=='sphere'):
            shapes.append([shape_args[i][0],shape_args[i][1],shape_args[i][2],shape_args[i][3]])
            print('Sphere arguments: inner_radius %.3f - outer_radius %.3f - num_radius %d' % \
              (shapes[i][1], shapes[i][2], shapes[i][3]))

        elif (shape_args[i][0]=='frustum') or (shape_args[i][0]=='cylinder'):
            if (shape_args[i][1][0]=='-'):
                flip = True
                if (shape_args[i][1][1:]=='minor'):
                    axis = 'disk minor axis'
                else:
                    axis = shape_args[i][1][1]
            elif (shape_args[i][1]=='minor'):
                flip = False
                axis = 'disk minor axis'
            else:
                flip = False
                axis = shape_args[i][1]

            if (shape_args[i][0]=='frustum'):
                shapes.append([shape_args[i][0], shape_args[i][2], shape_args[i][3], shape_args[i][4], axis, flip, shape_args[i][5]])
                if (flip):
                    print('Frustum arguments: axis - flipped %s inner_radius - %.3f outer_radius - %.3f num_steps - %d opening_angle - %d' % \
                      (axis, shapes[i][1], shapes[i][2], shapes[i][3], shapes[i][6]))
                else:
                    print('Frustum arguments: axis - %s inner_radius - %.3f outer_radius - %.3f num_steps - %d opening_angle - %d' % \
                      (str(axis), shapes[i][1], shapes[i][2], shapes[i][3], shapes[i][6]))

            if (shape_args[i][0]=='cylinder'):
                if (shape_args[i][5]!='height') and (shape_args[i][5]!='radius'):
                    sys.exit("I don't understand which way you want to calculate fluxes. Specify 'height' or 'radius'.")
                shapes.append([shape_args[i][0], shape_args[i][2], shape_args[i][3], shape_args[i][6], axis, flip, shape_args[i][4], shape_args[i][5]])
                if (flip):
                    print('Cylinder arguments: axis flipped %s - bottom_edge %.3f - top_edge %.3f - radius %.3f - step_direction %s - num_steps %d' % \
                      (axis, shapes[i][1], shapes[i][2], shapes[i][6], shapes[i][7], shapes[i][3]))
                else:
                    print('Cylinder arguments: axis %s - bottom_edge %.3f - top_edge %.3f - radius %.3f - step_direction %s - num_steps %d' % \
                      (str(axis), shapes[i][1], shapes[i][2], shapes[i][6], shapes[i][7], shapes[i][3]))

        else:
            sys.exit("That shape has not been implemented. Ask Cassi to add it.")

    # Check to make sure if anything is a cylinder, then there are no more shapes
    cyls = 0
    for i in range(len(shapes)):
        if (shapes[i][0]=='cylinder'):
            cyls += 1
    if (cyls > 0):
        if (cyls!=len(shapes)) or (cyls > 1):
            sys.exit("You can't have more than one cylinder or mix cylinders with other shapes! Calculate them separately.")

    # Check to make sure if multiple shapes are specified that they have the same inner_radius, outer_radius, and num_steps
    inner = shapes[0][1]
    outer = shapes[0][2]
    numsteps = shapes[0][3]
    for i in range(len(shapes)):
        if (shapes[i][1]!=inner) or (shapes[i][2]!=outer) or (shapes[i][3]!=numsteps):
            sys.exit('When specifying multiple shapes, you must give the same inner_radius,\n' + \
            'outer_radius, and num_steps for all shapes!')

    if (units_rvir):
        print('Shape arguments are in units of Rvir.')
    else:
        print('Shape arguments are in units of kpc.')

    return shapes

def segment_region(x, y, z, theta, phi, radius, shapes,
                   x_disk=False, y_disk=False, z_disk=False,
                   units_rvir=False, Rvir=100.):
    '''This function reads in arrays of x, y, z, theta_pos, phi_pos, and radius values and returns a
    boolean list of the same size that is True if a cell is contained within a shape in the list of
    shapes given by 'shapes' and is False otherwise. If disk-relative coordinates are needed for some
    shapes, they can be passed in with the optional x_disk, y_disk, z_disk.'''

    bool_inshape = np.zeros(len(x), dtype=bool)

    for i in range(len(shapes)):
        if (shapes[i][0]=='sphere'):
            if (units_rvir):
                inner_radius = shapes[i][1]*Rvir
                outer_radius = shapes[i][2]*Rvir
            else:
                inner_radius = shapes[i][1]
                outer_radius = shapes[i][2]
            bool_insphere = (radius > inner_radius) & (radius < outer_radius)
            bool_inshape = bool_inshape | bool_insphere
        elif (shapes[i][0]=='frustum'):
            if (units_rvir):
                inner_radius = shapes[i][1]*Rvir
                outer_radius = shapes[i][2]*Rvir
            else:
                inner_radius = shapes[i][1]
                outer_radius = shapes[i][2]
            op_angle = shapes[i][6]
            axis = shapes[i][4]
            flip = shapes[i][5]
            if (flip):
                min_theta = np.pi-op_angle*np.pi/180.
                max_theta = np.pi
            else:
                min_theta = 0.
                max_theta = op_angle*np.pi/180.
            if (axis=='x'):
                theta_frus = np.arccos(np.sin(theta)*np.cos(phi))
                phi_frus = np.arctan2(np.cos(theta), np.sin(theta)*np.sin(phi))
            if (axis=='y'):
                theta_frus = np.arccos(np.sin(theta)*np.sin(phi))
                phi_frus = np.arctan2(np.sin(theta)*np.cos(phi), np.cos(theta))
            if (axis=='disk minor axis'):
                theta_frus = np.arccos(z_disk/radius)
                phi_frus = np.arctan2(y_disk, x_disk)
            if (type(axis)==tuple) or (type(axis)==list):
                axis = np.array(axis)
                norm_axis = axis / np.sqrt((axis**2.).sum())
                # Define other unit vectors orthagonal to the angular momentum vector
                np.random.seed(99)
                x_axis = np.random.randn(3)            # take a random vector
                x_axis -= x_axis.dot(norm_axis) * norm_axis       # make it orthogonal to L
                x_axis /= np.linalg.norm(x_axis)            # normalize it
                y_axis = np.cross(norm_axis, x_axis)           # cross product with L
                x_vec = np.array(x_axis)
                y_vec = np.array(y_axis)
                z_vec = np.array(norm_axis)
                # Calculate the rotation matrix for converting from original coordinate system
                # into this new basis
                xhat = np.array([1,0,0])
                yhat = np.array([0,1,0])
                zhat = np.array([0,0,1])
                transArr0 = np.array([[xhat.dot(x_vec), xhat.dot(y_vec), xhat.dot(z_vec)],
                                     [yhat.dot(x_vec), yhat.dot(y_vec), yhat.dot(z_vec)],
                                     [zhat.dot(x_vec), zhat.dot(y_vec), zhat.dot(z_vec)]])
                rotationArr = np.linalg.inv(transArr0)
                x_rot = rotationArr[0][0]*np.sin(theta)*np.cos(phi) + rotationArr[0][1]*np.sin(theta)*np.sin(phi) + rotationArr[0][2]*np.cos(theta)
                y_rot = rotationArr[1][0]*np.sin(theta)*np.cos(phi) + rotationArr[1][1]*np.sin(theta)*np.sin(phi) + rotationArr[1][2]*np.cos(theta)
                z_rot = rotationArr[2][0]*np.sin(theta)*np.cos(phi) + rotationArr[2][1]*np.sin(theta)*np.sin(phi) + rotationArr[2][2]*np.cos(theta)
                theta_frus = np.arccos(z_rot)
                phi_frus = np.arctan2(y_rot, x_rot)
            bool_infrus = (theta_frus >= min_theta) & (theta_frus <= max_theta) & (radius >= inner_radius) & (radius <= outer_radius)
            bool_inshape = bool_inshape | bool_infrus
        elif (shapes[i][0]=='cylinder'):
            if (units_rvir):
                bottom_edge = shapes[i][1]*Rvir
                top_edge = shapes[i][2]*Rvir
                cyl_radius = shapes[i][6]*Rvir
            else:
                bottom_edge = shapes[i][1]
                top_edge = shapes[i][2]
                cyl_radius = shapes[i][6]                
            axis = shapes[i][4]
            flip = shapes[i][5]
            if (flip): mult = -1.
            else: mult = 1.
            if (axis=='z'):
                norm_coord = mult*z
                rad_coord = np.sqrt(x**2. + y**2.)
            if (axis=='x'):
                norm_coord = mult*x
                rad_coord = np.sqrt(y**2. + z**2.)
            if (axis=='y'):
                norm_coord = mult*y
                rad_coord = np.sqrt(x**2. + z**2.)
            if (axis=='disk minor axis'):
                norm_coord = mult*z_disk
                rad_coord = np.sqrt(x_disk**2. + y_disk**2.)
            if (type(axis)==tuple) or (type(axis)==list):
                axis = np.array(axis)
                norm_axis = axis / np.sqrt((axis**2.).sum())
                # Define other unit vectors orthagonal to the angular momentum vector
                np.random.seed(99)
                x_axis = np.random.randn(3)            # take a random vector
                x_axis -= x_axis.dot(norm_axis) * norm_axis       # make it orthogonal to L
                x_axis /= np.linalg.norm(x_axis)            # normalize it
                y_axis = np.cross(norm_axis, x_axis)           # cross product with L
                x_vec = np.array(x_axis)
                y_vec = np.array(y_axis)
                z_vec = np.array(norm_axis)
                # Calculate the rotation matrix for converting from original coordinate system
                # into this new basis
                xhat = np.array([1,0,0])
                yhat = np.array([0,1,0])
                zhat = np.array([0,0,1])
                transArr0 = np.array([[xhat.dot(x_vec), xhat.dot(y_vec), xhat.dot(z_vec)],
                                     [yhat.dot(x_vec), yhat.dot(y_vec), yhat.dot(z_vec)],
                                     [zhat.dot(x_vec), zhat.dot(y_vec), zhat.dot(z_vec)]])
                rotationArr = np.linalg.inv(transArr0)
                x_rot = rotationArr[0][0]*x + rotationArr[0][1]*y + rotationArr[0][2]*z
                y_rot = rotationArr[1][0]*x + rotationArr[1][1]*y + rotationArr[1][2]*z
                z_rot = rotationArr[2][0]*x + rotationArr[2][1]*y + rotationArr[2][2]*z
                norm_coord = mult*z_rot
                rad_coord = np.sqrt(x_rot**2. + y_rot**2.)
            bool_incyl = (norm_coord >= bottom_edge) & (norm_coord <= top_edge) & (rad_coord <= cyl_radius)
            bool_inshape = bool_inshape | bool_incyl
       
    if (shapes[0][0]=='cylinder') and (shapes[0][7]=='radius'):
        return bool_inshape, rad_coord
    elif (shapes[0][0]=='cylinder') and (shapes[0][7]=='height'):
        return bool_inshape, norm_coord
    else:
        return bool_inshape

