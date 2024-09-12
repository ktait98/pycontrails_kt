import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import dask.array as da
import pyproj
import sys
from joblib import Parallel, delayed
from tqdm import tqdm

pd.set_option("display.max_rows", 200)
import pdb

from pycontrails import Flight, Fleet, MetDataset
from pycontrails.core import models
from pycontrails.datalib.ecmwf import ERA5
from pycontrails.physics import geo, thermo, units, constants

# from pycontrails.models.ps_model import PSFlight
# from pycontrails.models.emissions import Emissions
from pycontrails.ext.flight_gen import FlightGen
from pycontrails.models.boxmodel.boxm_ac import Boxm

# from pycontrails.models.dry_advection import DryAdvection
from pycontrails.core.met_var import (
    AirTemperature,
    RelativeHumidity,
    SpecificHumidity,
    EastwardWind,
    NorthwardWind,
    VerticalVelocity,
)

# meteorological parameters
met_params = {
    "air_temperature": 235.0,  # K
    "specific_humidity": 0.003,  # 1
    "relative_humidity": 0.5,  # 1
    "eastward_wind": 0.0,  # m/s
    "northward_wind": 0.0,  # m/s
    "lagrangian_tendency_of_air_pressure": 0.0,  # m/s
}

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
    "max_age": pd.Timedelta(hours=2),  # maximum age of the plume
    "depth": 50.0,  # initial plume depth, [m]
    "width": 50.0,  # initial plume width, [m]
    "shear": 0.01,  # wind shear [1/s]
}

# chemistry sim parameters
chem_params = {
    "t0_chem": pd.to_datetime("2022-01-20 12:00:00"),  # chemistry start time
    "rt_chem": pd.Timedelta(hours=3),  # chemistry runtime
    "ts_chem": pd.Timedelta(seconds=20),  # chemistry time step
    "lat_bounds": (0.0, 1.0),  # lat bounds [deg]
    "lon_bounds": (0.0, 1.0),  # lon bounds [deg]
    "alt_bounds": (12000, 13000),  # alt bounds [m]
    "hres_pl": 0.01,  # horizontal resolution of the plume, [deg]
    "hres_chem": 0.01,  # horizontal resolution [deg]
    "vres_chem": 500,  # vertical resolution [m]
}

lats_pl = np.arange(
    chem_params["lat_bounds"][0], chem_params["lat_bounds"][1] + chem_params["hres_pl"], chem_params["hres_pl"]
)

lons_pl = np.arange(
    chem_params["lon_bounds"][0], chem_params["lon_bounds"][1] + chem_params["hres_pl"], chem_params["hres_pl"]
)

lats = np.arange(
    chem_params["lat_bounds"][0], chem_params["lat_bounds"][1] + chem_params["hres_chem"], chem_params["hres_chem"]
)

lons = np.arange(
    chem_params["lon_bounds"][0], chem_params["lon_bounds"][1] + chem_params["hres_chem"], chem_params["hres_chem"]
)

alts = np.arange(
    chem_params["alt_bounds"][0], chem_params["alt_bounds"][1] + chem_params["vres_chem"], chem_params["vres_chem"]
)

times = pd.date_range(
    start=chem_params["t0_chem"],
    end=chem_params["t0_chem"] + chem_params["rt_chem"],
    freq=chem_params["ts_chem"],
)

def boxm_run(lons, lons_pl, lats, lats_pl, alts, times, fl_params, plume_params, chem_params, met_params):
        
        met = gen_met(lons, lats, alts, times, met_params)

        fl_gen = FlightGen(met, fl_params, plume_params, chem_params)

        fl = fl_gen.traj_gen()

        fl = fl_gen.calc_fb_emissions()

        fl_df, pl_df = fl_gen.sim_plumes()

        emi = fl_gen.plume_to_grid(lats_pl, lons_pl, alts, times)

        boxm = Boxm(met=met, fl_params=fl_params, plume_params=plume_params, chem_params=chem_params)
        
        boxm.eval(emi)

        return fl_df, pl_df, boxm

def gen_met(lons, lats, alts, times, met_params):

    # generate artifical met dataset (boxm currently only supports zero-wind scenarios)
    data_vars = {
        param: (
            ["longitude", "latitude", "level", "time"],
            da.full(
                (len(lons), len(lats), len(alts), len(times)),
                value,
                chunks=(len(lons), len(lats), len(alts), 100),
            ),
        )
        for param, value in met_params.items()
    }

    met = xr.Dataset(
        data_vars,
        coords={"longitude": lons, "latitude": lats, "level": units.m_to_pl(alts), "time": times},
    )

    met = MetDataset(met)

    return met

# Run the box model
fl_df, pl_df, boxm = boxm_run(lons, lons_pl, lats, lats_pl, alts, times, fl_params, plume_params, chem_params, met_params)



# Plot the heatmap
fig1, ax1 = plt.subplots()
ax1.set_xticks(np.arange(chem_params["lon_bounds"][0], chem_params["lon_bounds"][1], 0.05))
ax1.set_yticks(np.arange(chem_params["lat_bounds"][0], chem_params["lat_bounds"][1], 0.05))


ts = 6

# Plot the heatmap
heatmap_data = boxm.boxm_ds_unstacked.emi.sel(emi_species="NO", time=fl_df.time[ts]).sel(level=178.6, method="nearest").transpose("latitude", "longitude")
heatmap_data.plot(ax=ax1, cmap='summer')  # You can choose a colormap of your preference


scat_fl = ax1.scatter(fl_df["longitude"].loc[fl_df["time"] == fl_df["time"][ts]],
                      fl_df["latitude"].loc[fl_df["time"] == fl_df["time"][ts]], 
                      s=5, c="red", label="Flight path")

scat_pl = ax1.scatter(pl_df["longitude"].loc[pl_df["time"] == fl_df["time"][ts]],
                      pl_df["latitude"].loc[pl_df["time"] == fl_df["time"][ts]],
                      s=10E-2 * pl_df["width"].loc[pl_df["time"] == fl_df["time"][ts]], c="blue", label="Plume evolution")

ax1.legend(loc="upper left")
ax1.set_xlim([chem_params["lon_bounds"][0], chem_params["lon_bounds"][1]])
ax1.set_ylim([chem_params["lat_bounds"][0], chem_params["lat_bounds"][1]])
#plt.grid()
plt.show()