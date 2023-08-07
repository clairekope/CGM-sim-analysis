#!/usr/bin/env python
# coding: utf-8

import yt
from yt.units import dimensions
from yt.fields.api import ValidateParameter # req yt 4 for correctness
import numpy as np

# Field Functions

def _overdensity(field, data):
    odthresh = data.get_field_parameter("overdensity_threshold")
    return (data[("gas","density")] > odthresh).astype(float)

def _velocity_divergence(field, data):
    return (data[("gas", "velocity_divergence")] < 0).astype(float)

def _freefall_time(field, data):
    total_density = data[("gas","density")] + data[("gas", "dark_matter_density")]
    return np.sqrt(3*np.pi / (32 * yt.units.gravitational_constant * total_density))

def _cold_or_cooling(field, data):
    tempthresh = data.get_field_parameter("temperature_threshold")
    cold = data[("gas", "temperature")] < tempthresh
    cooling = data[("gas", "cooling_time")] < data[("gas", "freefall_time")]
    return np.logical_or(cold, cooling).astype(float)

def _jeans_mass(field, data):
    """
    WARNING: Enzo uses a constant assumption for mean molecular weight,
    whereas this field does not.
    """
    return (data[("gas", "cell_mass")] > data[("gas","jeans_mass")]).astype(float)

# Add fields to yt

yt.add_field(
     ("index", "overdensity_criterion"),
     function = _overdensity,
     sampling_type = "cell",
     dimensions = dimensions.dimensionless,
     take_log = False,
     validators = [ValidateParameter(["overdensity_threshold"])]
)

yt.add_field(
    ("index", "velocity_divergence_criterion"),
    function = _velocity_divergence,
    sampling_type = "cell",
    dimensions = dimensions.dimensionless,
    take_log = False
)

yt.add_field(
    ("gas", "freefall_time"),
    function = _freefall_time,
    sampling_type = "cell",
    dimensions = dimensions.time,
    take_log = True,
)

yt.add_field(
    ("index", "temperature_criterion"),
    function = _cold_or_cooling,
    sampling_type = "cell",
    dimensions = dimensions.dimensionless,
    take_log = False,
    validators = [ValidateParameter(["temperature_threshold"])]
)

yt.add_field(
    ("index", "jeans_mass_criterion"),
    function = _jeans_mass,
    sampling_type = "cell",
    dimensions = dimensions.dimensionless,
    take_log = False
)

# Sample usage

if __name__ == "__main__":

    ds = yt.load("../original_sims/sample_data/fid/DD0030/DD0030")

    dsk = ds.disk('c', [0,0,1], (30,'kpc'), (3,'kpc'))

    dsk.set_field_parameter("overdensity_threshold", 1e2*yt.units.mh_cgs/yt.units.cm**3)
    dsk.set_field_parameter("temperature_threshold", 3e3*yt.units.K)
    print("Overdensity threshold:", dsk.field_parameters["overdensity_threshold"])
    print("Temperature threshold:", dsk.field_parameters["temperature_threshold"])

    slc = yt.SlicePlot(ds, 'z', center='c', 
             fields=[("index","overdensity_criterion"),
                     ("index","velocity_divergence_criterion"),
                     ("index","temperature_criterion"),
                     ("index","jeans_mass_criterion")],
             width=(30,'kpc'), data_source=dsk)

    slc.save()


