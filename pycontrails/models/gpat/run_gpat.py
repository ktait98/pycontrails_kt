#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import dask.array as da
import pyproj
import sys
from joblib import Parallel, delayed
from tqdm import tqdm
from pycontrails import Flight, Fleet, MetDataset
from pycontrails.core import models
from pycontrails.datalib.ecmwf import ERA5
from pycontrails.physics import geo, thermo, units, constants
from pycontrails.models.gpat.gpat import GPAT

# flight trajectory parameters
fl_params = {
    "t0_fl": pd.to_datetime("2022-01-20 13:00:00"),  # flight start time
    "rt_fl": pd.Timedelta(minutes=60),  # flight run time
    "ts_fl": pd.Timedelta(minutes=2),  # flight time step
    "ac_type": "A320",  # aircraft type
    "fl0_speed": 100.0,  # m/s
    "fl0_heading": 0.0,  # deg
    "fl0_coords0": (0.1, 0.125, 12500),  # lat, lon, alt [deg, deg, m]
    "sep_dist": (5000, 2000, 0),  # dx, dy, dz [m]
    "n_ac": 1,  # number of aircraft
}

# plume dispersion parameters
plume_params = {
    "dt_integration": pd.Timedelta(minutes=2),  # integration time step
    "max_age": pd.Timedelta(hours=1),  # maximum age of the plume
    "depth": 50.0,  # initial plume depth, [m]
    "width": 50.0,  # initial plume width, [m]
    "shear": 0.01,  # wind shear [1/s]
    "hres_pl": 0.02, # horizontal resolution of the plume [deg]
    "vres_pl": 500 # vertical resolution of the plume [m]
}

# chemistry sim parameters
sim_params = {
    "t0_sim": pd.to_datetime("2022-01-20 12:00:00"),  # chemistry start time
    "rt_sim": pd.Timedelta(days=5),  # chemistry runtime
    "ts_sim": pd.Timedelta(seconds=20),  # chemistry time step
    "lat_bounds": (0.0, 1.0),  # lat bounds [deg]
    "lon_bounds": (0.0, 1.0),  # lon bounds [deg]
    "alt_bounds": (12000, 13000),  # alt bounds [m]
    "hres_sim": 0.02,  # horizontal resolution [deg]
    "vres_sim": 500,  # vertical resolution [m]
    "eastward_wind": 0.0,  # m/s
    "northward_wind": 0.0,  # m/s
    "lagrangian_tendency_of_air_pressure": 0.0,  # m/s
    "species_out": np.array([3, 4, 5, 6, 8, 9]) #, 14, 21, 22, 39, 198])
}

gpat = GPAT(fl_params, plume_params, sim_params)

gpat.eval()