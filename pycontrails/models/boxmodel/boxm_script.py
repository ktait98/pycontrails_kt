import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

pd.set_option("display.max_rows", 200)

from pycontrails import MetDataset

# from pycontrails.models.ps_model import PSFlight
# from pycontrails.models.emissions import Emissions
from pycontrails.ext.flight_gen import FlightGen
from pycontrails.models.boxmodel.boxm import Boxm
from pycontrails.physics import units

# from pycontrails.models.dry_advection import DryAdvection

# meteorological parameters
met_params = {
    "air_temperature": 240.0,  # K
    "specific_humidity": 0.001,  # 1
    "relative_humidity": 0.5,  # 1
    "eastward_wind": 0.0,  # m/s
    "northward_wind": 0.0,  # m/s
    "lagrangian_tendency_of_air_pressure": 0.0,  # m/s
}

# flight trajectory parameters
fl_params = {
    "t0_fl": pd.to_datetime("2022-03-02 08:00:00"),  # flight start time
    "rt_fl": pd.Timedelta(hours=1),  # flight run time
    "ts_fl": pd.Timedelta(minutes=2),  # flight time step
    "ac_type": "A320",  # aircraft type
    "fl0_speed": 100.0,  # m/s
    "fl0_heading": 45.0,  # deg
    "fl0_coords0": (-0.8, -0.8, 11500),  # lat, lon, alt [deg, deg, m]
    "sep_dist": (5000, 2000, 0),  # dx, dy, dz [m]
    "n_ac": 2,  # number of aircraft
}

# plume dispersion parameters
plume_params = {
    "dt_integration": pd.Timedelta(minutes=5),  # integration time step
    "max_age": pd.Timedelta(hours=1),  # maximum age of the plume
    "depth": 50.0,  # initial plume depth, [m]
    "width": 40.0,  # initial plume width, [m]
    "shear": 0.01,  # wind shear [1/s]
}

# chemistry sim parameters
chem_params = {
    "t0_chem": pd.to_datetime("2022-03-02 07:00:00"),  # chemistry start time
    "rt_chem": pd.Timedelta(days=1),  # chemistry runtime
    "ts_chem": pd.Timedelta(minutes=5),  # chemistry time step
    "lat_bounds": (-1.0, 0.0),  # lat bounds [deg]
    "lon_bounds": (-1.0, 0.0),  # lon bounds [deg]
    "alt_bounds": (11000, 12500),  # alt bounds [m]
    "hres_chem": 0.2,  # horizontal resolution [deg]
    "vres_chem": 500,  # vertical resolution [m]
}

# create lists for lats, lons, alts, and times based on chem params
lats = np.arange(
    chem_params["lat_bounds"][0], chem_params["lat_bounds"][1], chem_params["hres_chem"]
)

lons = np.arange(
    chem_params["lon_bounds"][0], chem_params["lon_bounds"][1], chem_params["hres_chem"]
)

alts = np.arange(
    chem_params["alt_bounds"][0], chem_params["alt_bounds"][1], chem_params["vres_chem"]
)

times = pd.date_range(
    start=chem_params["t0_chem"],
    end=chem_params["t0_chem"] + chem_params["rt_chem"],
    freq=chem_params["ts_chem"],
)

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

# generate flight trajectories
fl_gen = FlightGen(met, fl_params, plume_params, chem_params)

flights = fl_gen.traj_gen()

# visualise the fleet
ax = plt.axes()
ax.set_xlim([lons[0], lons[-1]])
ax.set_ylim([lats[0], lats[-1]])
for fl in flights:
    fl.plot(ax=ax)

# estimate fuel burn and emissions using ps_model and emissions model
flights = fl_gen.calc_fb_emissions()


# simulate plume dispersion/advection using dry advection model
fl_df, pl_df = fl_gen.sim_plumes()

print(pl_df["time"].values)

# visualise plumes
fl_gen.anim_fl(fl_df, pl_df)

# convert plume dataframe to EMI geospatial xarray dataset
emi = fl_gen.plume_to_grid()

# init boxm simulation and generate chemistry dataset
boxm = Boxm(met=met, params=chem_params)

# run boxm simulation
chem = boxm.eval(emi)
chem.to_netcdf("chem.nc")


