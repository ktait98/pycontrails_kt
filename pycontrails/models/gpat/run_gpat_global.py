#!/usr/bin/env python
import numpy as np
import pandas as pd
import xarray as xr
from pycontrails.models.gpat.gpat_clean import GPAT_clean, SimParams, dict_to_dataclass
from dataclasses import asdict
from pycontrails.core import MetDataset
from pycontrails.physics import geo, thermo, units, constants
import pickle
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys

if len(sys.argv) != 2:
    print("Usage: run_gpat_global.py <month>")
    sys.exit(1)

month = sys.argv[1]

sim_params = {
    "t0_sim": pd.to_datetime("2022-01-01 12:00:00"),  # chemistry start time
    "rt_sim": pd.Timedelta(days=5),  # chemistry runtime
    "ts_sim": pd.Timedelta(seconds=20),  # chemistry time step
    "lat_bounds": (-87.5, 87.5),  # lat bounds [deg]
    "lon_bounds": (-177.5, 177.5),  # lon bounds [deg]
    "alt_bounds": (8000, 13000),  # alt bounds [m]
    "hres_sim": 5,  # horizontal resolution [deg]
    "vres_sim": 1000,  # vertical resolution [m]
    "eastward_wind": 0.0,  # m/s
    "northward_wind": 0.0,  # m/s
    "lagrangian_tendency_of_air_pressure": 0.0,  # m/s
    "species_in": ("NO", "NO2", "CO", "HCHO", "CH3CHO", "C2H4", "C3H6", "C2H2", "BENZENE"),
    "species_out": ("O3", "NO2", "NO", "NO3", "HNO3", "PAN", "HONO", "HO2", "OH","H2O2", 
                    "CO", "HCHO", "CH4"),
    "gpat_path": os.getcwd() +"/",
    "job_id": f"global_{month}",
    "run_gpat": None
}

sim_params = dict_to_dataclass(SimParams, sim_params)

print("SimParams:", asdict(sim_params)) 

gpat = GPAT_clean(sim_params)

gpat.met = gpat.gen_met()
gpat.bg_chem = gpat.gen_bg_chem()

gpat.chem = gpat.run_boxm()

gpat.gen_outputs()