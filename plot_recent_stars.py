#!/usr/bin/env python
# coding: utf-8
import yt
import numpy as np

@yt.particle_filter(requires=["particle_type"], filtered_type='all')
def recent_stars(pfilter, data):
    
    dt = data.ds.quan(data.ds.index.parameters['dtDataDump'], 'code_time')
    min_creation_time = data.ds.current_time - dt
    
    filter1 = data[(pfilter.filtered_type, "particle_type")] == 2
    filter2 = data[(pfilter.filtered_type, "creation_time")] > min_creation_time
    
    filter = np.logical_and(filter1, filter2)
    
    return filter

datasets = yt.load("DD????/DD????")
for ds in datasets.piter():

    ds.add_particle_filter('recent_stars')

    center = ds.quan(0.5,'code_length')
    rs = ds.quan(3.5,'kpc')

    width = ds.quan(100,'kpc')
    thickness = rs

    # Edge on
    rect = ds.region([center, center, center],
                     [center-thickness/2, center-width/2, center-width/2],
                     [center+thickness/2, center+width/2, center+width/2])
    assert (rect.get_field_parameter("center") == center).all()

    p1 = yt.ProjectionPlot(ds, 'x', 'density', width=width,
                           data_source=rect, weight_field ='ones')
    p1.set_zlim("density", 1e-32, 1e-24)
    p1.annotate_particles(width=width, ptype='recent_stars', p_size=5)
    p1.save(ds.basename+"_edge_on_stars.png")

    # Face on
    rect = ds.region([center, center, center],
                     [center-width/2, center-width/2, center-thickness/2],
                     [center+width/2, center+width/2, center+thickness/2])
    assert (rect.get_field_parameter("center") == center).all()

    p2 = yt.ProjectionPlot(ds, 'z', 'density', width=width,
                           data_source=rect, weight_field ='ones')
    p2.set_zlim("density", 1e-29, 1e-23)
    p2.annotate_particles(width=width, ptype='recent_stars', p_size=5)
    p2.save(ds.basename+"_face_on_stars.png")

