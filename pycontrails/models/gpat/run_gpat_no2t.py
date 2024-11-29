#!/usr/bin/env python

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from pycontrails.models.gpat.gpat import GPAT, FlParams, PlumeParams, SimParams, parse_args, update_fl_params_from_args, update_plume_params_from_args, update_sim_params_from_args, dict_to_dataclass
from dataclasses import asdict
import os

fl_params = {
    "t0_fl": pd.to_datetime("2022-01-20 13:00:00"),  # flight start time
    "rt_fl": pd.Timedelta(minutes=60),  # flight run time
    "ts_fl": pd.Timedelta(minutes=2),  # flight time step
    "ac_type": "A320",  # aircraft type
    "fl0_speed": 100.0,  # m/s
    "fl0_heading": 45.0,  # deg
    "fl0_coords0": (47.1, -32.9, 12500),  # lat, lon, alt [deg, deg, m]
    "sep_dist": (0, 0, 0),  # dx, dy, dz [m]
    "n_ac": 0,  # number of aircraft
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


fl_params = dict_to_dataclass(FlParams, fl_params)
plume_params = dict_to_dataclass(PlumeParams, plume_params)

locations = [
    {
        "name": "NA",
        "lat": 47.5,
        "lon": -32.5,
        "alt": 12500,
        "datetime": "20-01-2022 12:00:00"
    },
    {
        "name": "US",
        "lat": 37.5,
        "lon": -97.5,
        "alt": 11500,
        "datetime": "10-11-2022 12:00:00"
    },
    {
        "name": "EU",
        "lat": 42.5,
        "lon": 7.5,
        "alt": 9500,
        "datetime": "20-05-2023 12:00:00"
    },
    {
        "name": "SEA",
        "lat": 22.5,
        "lon": 102.5,
        "alt": 10500,
        "datetime": "05-03-2022 12:00:00"
    },
    {
        "name": "SA",
        "lat": -27.5,
        "lon": -67.5,
        "alt": 13500,
        "datetime": "15-08-2022 12:00:00"

    }
]

fig, ax = plt.subplots()

# Example of looping through the data
for location in locations:
    location['NO2_equiv'] = []
    location['del_NO2'] = []

    for month in range(1,13):
        print(f"Month: {month}")

        if month < 10:

            sim_params = {
                "t0_sim": pd.to_datetime(f"2022-0{month}-20 12:00:00"),  # chemistry start time
                "rt_sim": pd.Timedelta(days=5),  # chemistry runtime
                "ts_sim": pd.Timedelta(seconds=20),  # chemistry time step
                "lat_bounds": (location['lat'] - 0.5, location['lat'] + 0.5),  # lat bounds [deg]
                "lon_bounds": (location['lon'] - 0.5, location['lon'] + 0.5),  # lon bounds [deg]
                "alt_bounds": (location['alt'] - 500, location['alt'] + 500),  # alt bounds [m]
                "hres_sim": 0.5,  # horizontal resolution [deg]
                "vres_sim": 500,  # vertical resolution [m]
                "eastward_wind": 0.0,  # m/s
                "northward_wind": 0.0,  # m/s
                "lagrangian_tendency_of_air_pressure": 0.0,  # m/s
                "species_in": ("NO", "NO2", "CO", "HCHO", "CH3CHO", "C2H4", "C3H6", "C2H2", "BENZENE"),
                "species_out": ("O3", "NO2", "NO", "NO3", "HNO3", "PAN", "HONO", "HO2", "OH","H2O2", 
                                "CO", "HCHO", "CH4"),
                "gpat_path": os.getcwd() +"/",
                "job_id": None,
                "run_gpat": None
            }
        else:
            sim_params = {
                "t0_sim": pd.to_datetime(f"2022-{month}-20 12:00:00"),  # chemistry start time
                "rt_sim": pd.Timedelta(days=5),  # chemistry runtime
                "ts_sim": pd.Timedelta(seconds=20),  # chemistry time step
                "lat_bounds": (location['lat'] - 0.5, location['lat'] + 0.5),  # lat bounds [deg]
                "lon_bounds": (location['lon'] - 0.5, location['lon'] + 0.5),  # lon bounds [deg]
                "alt_bounds": (location['alt'] - 500, location['alt'] + 500),  # alt bounds [m]
                "hres_sim": 0.5,  # horizontal resolution [deg]
                "vres_sim": 500,  # vertical resolution [m]
                "eastward_wind": 0.0,  # m/s
                "northward_wind": 0.0,  # m/s
                "lagrangian_tendency_of_air_pressure": 0.0,  # m/s
                "species_in": ("NO", "NO2", "CO", "HCHO", "CH3CHO", "C2H4", "C3H6", "C2H2", "BENZENE"),
                "species_out": ("O3", "NO2", "NO", "NO3", "HNO3", "PAN", "HONO", "HO2", "OH","H2O2", 
                                "CO", "HCHO", "CH4"),
                "gpat_path": os.getcwd() +"/",
                "job_id": None,
                "run_gpat": None
            }


        sim_params = dict_to_dataclass(SimParams, sim_params)


        gpat = GPAT(fl_params, plume_params, sim_params)

        gpat.met = gpat.gen_met()
        gpat.bg_chem = gpat.gen_bg_chem()

        # gpat.met.data['latitude'] = gpat.met.data['latitude'] + 2.5
        # gpat.met.data['longitude'] = gpat.met.data['longitude'] + 2.5

        # gpat.bg_chem['latitude'] = gpat.bg_chem['latitude'] + 2.5
        # gpat.bg_chem['longitude'] = gpat.bg_chem['longitude'] + 2.5

        NO = gpat.bg_chem["bg_chem"].sel(species="NO")
        NO2 = gpat.bg_chem["bg_chem"].sel(species="NO2")
        CO = gpat.bg_chem["bg_chem"].sel(species="CO")
        CH4 = gpat.bg_chem["bg_chem"].sel(species="CH4")
        OH = gpat.bg_chem["bg_chem"].sel(species="OH")
        T = gpat.met["air_temperature"].data.isel(time=0)
        M = gpat.met["M"].data.isel(time=0)

        # Rate coefficients
        k1 = 5.4e-14 * (T / 298)**1.5 * np.exp(250 / T)
        k7 = 2.45e-12 * np.exp(-1775 / T)
        k0 = 2.5e-30 * (T / 300)**-4.4
        k_alpha = 1.6e-11 * (T / 300)**-1.7

        # Calculate k6[M]
        k6_M = (k0 * M) / (1 + (k0 * M) / k_alpha) * 0.6**(1 / (1 + np.log10(k0 * M / k_alpha)**2))

        NO2_equiv = (k1 * CO + k7 * CH4) / k6_M
        del_NO2 = NO2_equiv - NO2

        # Calculate [NO2]_equiv
        location['NO2_equiv'].append(NO2_equiv.sel(latitude=location['lat'], longitude=location['lon']).isel(level=2).item())
        location['del_NO2'].append(del_NO2.sel(latitude=location['lat'], longitude=location['lon']).isel(level=2).item())
    print(f"NO2_equiv: {location['NO2_equiv']}")
    print(f"del NO2: {location['del_NO2']}")
    location['del_NO2'] = xr.DataArray(location['del_NO2'], dims="time", coords={"time": pd.date_range(start="2022-01-01", periods=len(location['del_NO2']), freq="M")})
    # plot time series for each location on the same plot
    location['del_NO2'].plot(ax=ax, label=location['name'])

# Add labels and legend
ax.set_xlabel('Month')
ax.set_ylabel('Concentration [ppbv]')
ax.set_title('del NO2')
ax.legend()

plt.show()

