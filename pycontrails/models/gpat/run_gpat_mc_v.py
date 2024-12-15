#!/usr/bin/env python

import numpy as np
import pandas as pd
from pycontrails.models.gpat.gpat import GPAT, FlParams, PlumeParams, SimParams, parse_args, update_fl_params_from_args, update_plume_params_from_args, update_sim_params_from_args, dict_to_dataclass
from dataclasses import asdict
import os

# flight trajectory parameters
fl_params = {
    "t0_fl": pd.to_datetime("2022-01-20 13:00:00"),  # flight start time
    "rt_fl": pd.Timedelta(minutes=60),  # flight run time
    "ts_fl": pd.Timedelta(minutes=2),  # flight time step
    "ac_type": "A320",  # aircraft type
    "fl0_speed": 100.0,  # m/s
    "fl0_heading": 90.0,  # deg
    "fl0_coords0": (0.5, 0.1, 10500),  # lat, lon, alt [deg, deg, m]
    "sep_dist": (1000, 0, 0),  # dx, dy, dz [m]
    "n_ac": 2,  # number of aircraft
}

# plume dispersion parameters
plume_params = {
    "dt_integration": pd.Timedelta(minutes=2),  # integration time step
    "max_age": pd.Timedelta(hours=2),  # maximum age of the plume
    "depth": 50.0,  # initial plume depth, [m]
    "width": 50.0,  # initial plume width, [m]
    "shear": 0.01,  # wind shear [1/s]
    "hres_pl": 0.05, # horizontal resolution of the plume [deg]
    "vres_pl": 500, # vertical resolution of the plume [m]
    "n_slices": 10,  # number of plume slices
}

# chemistry sim parameters
sim_params = {
    "t0_sim": pd.to_datetime("2022-01-20 12:00:00"),  # chemistry start time
    "rt_sim": pd.Timedelta(hours=12),  # chemistry runtime
    "ts_sim": pd.Timedelta(seconds=20),  # chemistry time step
    "lat_bounds": (0.0, 1.0),  # lat bounds [deg]
    "lon_bounds": (0.0, 1.0),  # lon bounds [deg]
    "alt_bounds": (10000, 11000),  # alt bounds [m]
    "hres_sim": 0.05,  # horizontal resolution [deg]
    "vres_sim": 500,  # vertical resolution [m]
    "eastward_wind": 0.0,  # m/s
    "northward_wind": 0.0,  # m/s
    "lagrangian_tendency_of_air_pressure": 0.0,  # m/s
    "species_in": ("NO", "NO2", "CO", "HCHO", "CH3CHO", "C2H4", "C3H6", "C2H2", "BENZENE"),
    "species_out": ("O3", "NO2", "NO", "NO3", "HNO3", "PAN", "HONO", "HO2", "OH","H2O2", 
                    "CO", "HCHO", "CH4", "CH3O2"),
    "run_path": os.getcwd() +"/",
    "data_path": os.getcwd() +"/", # "/projects/Impact_of_aviation_on_climate
    "job_id": None,
    "run_gpat": None
}

fl_params = dict_to_dataclass(FlParams, fl_params)
plume_params = dict_to_dataclass(PlumeParams, plume_params)
sim_params = dict_to_dataclass(SimParams, sim_params)

updated_args = parse_args()

update_fl_params_from_args(fl_params, updated_args)
print("FlParams:", asdict(fl_params))

update_plume_params_from_args(plume_params, updated_args)
print("PlumeParams:", asdict(plume_params))

update_sim_params_from_args(sim_params, updated_args)

print("SimParams:", asdict(sim_params)) 

gpat = GPAT(fl_params, plume_params, sim_params)

if gpat.sim_params.run_gpat:
    gpat.eval()
else:
    print("GPAT simulation is not run.")
    print(f'Job ID is : {gpat.sim_params.job_id}')
