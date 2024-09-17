"""A photochemical trajectory model for simulating atmospheric chemistry and NOx emissions."""
from __future__ import annotations

import os
import pathlib
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Any
import subprocess

from pycontrails.physics import units
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import dask.array as da
from matplotlib.animation import FuncAnimation, PillowWriter

from pycontrails.core.met import MetDataset
from pycontrails.core.met_var import (
    AirPressure,
    AirTemperature,
    EastwardWind,
    NorthwardWind,
    RelativeHumidity,
    SpecificHumidity,
)

from pycontrails.core.models import Model, ModelParams
from pycontrails.models.boxmodel import boxm
from pycontrails.physics import constants, geo

### Boxm Model Parameters ###
@dataclass
class BoxmParams(ModelParams):
    """Default box model parameters."""

    lat_bounds: tuple[float, float] | None = None  # latmin, latmax [deg]
    lon_bounds: tuple[float, float] | None = None  # lonmin, lonmax [deg]
    alt_bounds: tuple[float, float] | None = None  # altmin, altmax [m]
    hres_pl: float = 0.05  # plume horizontal resolution [deg]
    hres_chem: float = 0.5  # chemistry horizontal resolution [deg]
    vres_chem: float = 500  # chemistry vertical resolution [m]

    t0_chem: datetime | None = None  # chemistry start time
    rt_chem: timedelta | None = None  # chemistry runtime
    ts_chem: timedelta | None = None  # chemistry timestep

### Boxm Model Class ###
class Boxm(Model):
    """
    Box model class definition.

    Compute evolution of chemical concentrations using boxm based on
    tropospheric chemistry scheme STOCHEM-CRI v2.2 R5 [ref].
    """

    name = "boxm"
    long_name = "Photochemical Trajectory Model"
    met_variables = (
        AirTemperature,
        EastwardWind,
        NorthwardWind,
        SpecificHumidity,
        RelativeHumidity,
        AirPressure,
    )
    default_params = BoxmParams
    met: MetDataset
    bg_chem: xr.Dataset
    met_required = True

    def __init__(
        self,
        met: MetDataset,
        fl_params: dict, plume_params: dict, chem_params: dict,
        **params_kwargs: Any,
    ) -> None:

        # call Model init
        super().__init__(met, params=chem_params, **params_kwargs)

        # check if met and chem datasets are provided
        if met:
            self.met = met

        self.bg_chem = grab_bg_chem(met)

        self.fl_params = fl_params
        self.plume_params = plume_params
        self.chem_params = chem_params
        
        # need to add tests to see if all variables are there and if they are
        # within the bounds of the model

    def eval(
        self,
        emi: MetDataset,
        **params: Any,
    ) -> MetDataset:
        """
        Eval function for the box model.

        Evaluate the photochemical trajectory model to compute the evolution
        of chemical concentrations, subject to emissions.
        """
        self.path = "/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/"
        self.file_path = "/home/ktait98/pycontrails_kt/pycontrails/models/files/"
        self.update_params(params)
        self.set_source(emi)
        self.source = self.require_source_type(MetDataset)

        self.process_datasets()
        self.init_boxm_ds()
        self.to_netcdfs()
        self.run_boxm()
        self.unstack()

    def process_datasets(self):
        met = self.met

        """Process the met, bg_chem, and emi datasets to prepare them for the box model."""

        # chunk met and emi data
        self.met.data = self.met.data.chunk(
            {"longitude": "auto", "latitude": "auto", "level": "auto", "time": "auto"}
        )

        self.source.data = self.source.data.chunk(
            {"longitude": "auto", "latitude": "auto", "level": "auto", "time": "auto"}
        )

        met = calc_M_H2O(met)

        met.data["sza"] = (
            ("latitude", "longitude", "time"),
            calc_sza(
                met["latitude"].data.values, met["longitude"].data.values, met["time"].data.values
            ),
        )

        # interpolate and downselect bg_chem to chem grid
        self.bg_chem = self.bg_chem.interp(
            longitude=met["longitude"].data,
            latitude=met["latitude"].data,
            level=met["level"].data,
            method="linear",
        ) 

        # interpolate and downselect emi to chem grid
        self.source.data = self.source.data.interp(
            longitude=met["longitude"].data,
            latitude=met["latitude"].data,
            level=met["level"].data,
            method="linear",
        )

    def init_boxm_ds(self):
        self.boxm_ds = xr.merge([self.met.data, self.bg_chem, self.source.data])

        self.boxm_ds = self.boxm_ds.drop_vars(
            [
                "specific_humidity",
                "relative_humidity",
                "eastward_wind",
                "northward_wind",
                "lagrangian_tendency_of_air_pressure",
                "month",
            ]
        )

        self.boxm_ds = self.boxm_ds.assign_attrs(dts=self.params["ts_chem"].total_seconds())

        self.boxm_ds["J"] = (["time", "level", "longitude", "latitude", "photol_params"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], 5)))

        self.boxm_ds["DJ"] = (["time", "level", "longitude", "latitude", "photol_coeffs"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], 5)))

        self.boxm_ds["RC"] = (["time", "level", "longitude", "latitude", "therm_coeffs"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], 5)))

        self.boxm_ds["Y"] = (["time", "level", "longitude", "latitude", "species"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], self.boxm_ds.dims["species"])))

    def to_netcdfs(self):
        """Convert the met, bg_chem, and emi datasets to boxm_ds.nc for use in the box model."""

        # stack datasets to get cell index for fortran
        self.boxm_ds_stacked = self.boxm_ds.stack(
            {"cell": ["level", "longitude", "latitude"]}
        )
        print(self.boxm_ds_stacked.indexes["cell"])

        self.boxm_ds_stacked = self.boxm_ds_stacked.reset_index("cell")

        

        # Delete any existing netCDF files
        if pathlib.Path(self.file_path + "boxm_ds.nc").exists():
            print("deleting boxm_ds.nc")
            pathlib.Path(self.file_path + "boxm_ds.nc").unlink()

        if pathlib.Path(self.file_path + "boxm_input.txt").exists():
            print("deleting boxm_input.txt")
            pathlib.Path(self.file_path + "boxm_input.txt").unlink()

            # delete any input and output files
        if pathlib.Path(self.file_path + "BACKITNE.OUT").exists():
            print("deleting BACKITNE.OUT")
            pathlib.Path(self.file_path + "BACKITNE.OUT").unlink()

        if pathlib.Path(self.file_path + "ZEN.OUT").exists():
            print("deleting ZEN.OUT")
            pathlib.Path(self.file_path + "ZEN.OUT").unlink()

        if pathlib.Path(self.file_path + "J.OUT").exists():
            print("deleting J.OUT")
            pathlib.Path(self.file_path + "J.OUT").unlink()

        if pathlib.Path(self.file_path + "DJ.OUT").exists():
            print("deleting DJ.OUT")
            pathlib.Path(self.file_path + "DJ.OUT").unlink()

        if pathlib.Path(self.file_path + "RC.OUT").exists():
            print("deleting RC.OUT")
            pathlib.Path(self.file_path + "RC.OUT").unlink()

        if pathlib.Path(self.file_path + "Y.OUT").exists():
            print("deleting Y.OUT")
            pathlib.Path(self.file_path + "Y.OUT").unlink()

        if pathlib.Path(self.file_path + "FL.OUT").exists():
            print("deleting FL.OUT")
            pathlib.Path(self.file_path + "FL.OUT").unlink()

        # Convert DataFrames to Datasets and write to netCDF
        self.boxm_ds_stacked.to_netcdf(self.file_path + "boxm_ds.nc", mode="w")

    def run_boxm(self):
        """Run the box model in fortran using subprocess."""
        
        subprocess.call(
            [self.path + "boxm"]
        )
        
        # open nc file
        self.boxm_ds = xr.open_dataset(self.file_path + "boxm_ds.nc")
        print(self.boxm_ds)

    def unstack(self):
        """Unstack the box model dataset."""
        
        # Convert the dataset to a Dask dataset
        self.boxm_ds = self.boxm_ds.chunk({'cell': 100})  # Adjust chunk size based on your memory
        
        # Convert 'level', 'lat', and 'lon' to coordinates
        self.boxm_ds_unstacked = self.boxm_ds.set_coords(['level', 'longitude', 'latitude'])
        print("coords set")
        
        # Create a multi-index for the 'cell' dimension
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.set_index(cell=['level', 'longitude', 'latitude'])
        
        # Unstack the dataset
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.unstack("cell")
        
        # # Compute the result to trigger the lazy evaluation
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.compute()

    def anim_chem(self, var1, var2, level, resample_freq='2min'):
        """Animate the chemical concentrations."""
        fig, (ax, cbar_ax) = plt.subplots(
            1, 2, gridspec_kw={"width_ratios": (0.9, 0.05), "wspace": 0.2}, figsize=(12, 8)
        )

        if var1 == "Y":
            boxm_da = self.boxm_ds_unstacked[var1].sel(species=var2).sel(level=level, method="nearest")

        if var1 == "emi":
            boxm_da = self.boxm_ds_unstacked[var1].sel(emi_species=var2).sel(level=level, method="nearest")

        if var1 == "J":
            boxm_da = self.boxm_ds_unstacked[var1].sel(photol_params=var2).sel(level=level, method="nearest")

        if var1 == "DJ":
            boxm_da = self.boxm_ds_unstacked[var1].sel(photol_coeffs=var2).sel(level=level, method="nearest")

        if var1 == "RC":
            boxm_da = self.boxm_ds_unstacked[var1].sel(therm_coeffs=var2).sel(level=level, method="nearest")

        print(boxm_da)

        times = boxm_da["time"].values
        times_resampled = pd.to_datetime(times).to_series().resample(resample_freq).asfreq().dropna().index

        print(f"New number of frames: {len(times_resampled)}")
       
        def heatmap_func(t):
            ax.cla()
            ax.set_title(t)

            boxm_da.sel(time=t).transpose("latitude", "longitude").plot(
                ax=ax, cbar_kwargs={"cax": cbar_ax}, add_colorbar=True, vmin=boxm_da.min(), vmax=boxm_da.max()
            )

        anim = FuncAnimation(fig, heatmap_func, frames=times_resampled, blit=False)

        filename = pathlib.Path(self.file_path + var1 + "_" + var2 + ".gif")

        anim.save(filename, dpi=300, writer=PillowWriter(fps=8))

    def mc_test(self, fl_df):
        """Check if mass is conserved in the box model."""

        mm = [30.01, 46.01, 28.01, 30.03, 44.05, 28.05, 42.08, 26.04, 78.11]  # g/mol
        NA = 6.022e23  # Avogadro's number

        total_vector_mass = 0
        total_grid_mass = 0
        for ts, time in enumerate(fl_df["time"][:-1]):

            # grab vector data
            for s, emi_species in enumerate(["NO"]):# self.boxm_ds_unstacked["emi_species"].data):
                
                vector_mass = fl_df[emi_species][ts] 
                # \
                #         * self.plume_params["width"] \
                #         * self.chem_params["vres_chem"] \
                #         * fl_df["true_airspeed"][ts] \
                #         * self.plume_params["dt_integration"].seconds
            
                total_vector_mass += vector_mass

                # grab plume mass from grid data
                grid_concs = self.boxm_ds_unstacked["emi"].sel(emi_species=emi_species, time=time).sel(level=178.6, method="nearest")

                grid_concs_over_zero = grid_concs.where(grid_concs > 0, drop=True)

                grid_mass = grid_concs_over_zero \
                    * self.boxm_ds_unstacked["M"].sel(time=time).sel(level=178.6, method="nearest") \
                    * 1e-9 \
                    * (mm[s] / NA) \
                    * self.chem_params["vres_chem"] \
                    * units.latitude_distance_to_m(self.chem_params["hres_chem"]) \
                    * units.longitude_distance_to_m(self.chem_params["hres_chem"], (self.chem_params["lat_bounds"][0] + self.chem_params["lat_bounds"][1]) / 2) \
                    * 1E+03 # convert to kg/m^3
                    

                grid_mass_sum = grid_mass.sum().values

                total_grid_mass += grid_mass_sum
            
            print(total_vector_mass, total_grid_mass)
            
            

        


    #     for emi_species in self.boxm_ds_unstacked["emi_species"].data:
    #         emi_m3_fl = fl_df[emi_species]


    #         emi_m3_boxm = self.boxm_ds_unstacked["emi"].sel(emi_species=emi_species).sum()

    #         print(f"Mass emitted in {emi_species}: {emi} g")
    #         print(f"Mass emitted in {emi_species}: {emi_m3} g")
        

### functions used in boxm ###

# calculate solar zenith angle
def calc_sza(latitudes, longitudes, timesteps):
    """Calculate szas for each cell at all timesteps."""
    sza = np.zeros((len(latitudes), len(longitudes), len(timesteps)))

    for lon, lonval in enumerate(longitudes):
        for lat, latval in enumerate(latitudes):

            theta_rad = geo.orbital_position(timesteps)

            sza[lat, lon, :] = np.arccos(
                geo.cosine_solar_zenith_angle(lonval, latval, timesteps, theta_rad)
            )
    return sza

# calculate number density of air molecules and H2O
def calc_M_H2O(met):
    """Calculate number density of air molecules at each pressure level M."""
    N_A = 6.022e23  # Avogadro's number

    # Get air density from pycontrails physics.thermo script
    rho_d = met["air_pressure"].data / (constants.R_d * met["air_temperature"].data)

    # Calculate number density of air (M) to feed into box model calcs
    met.data["M"] = (N_A / constants.M_d) * rho_d * 1e-6  # [molecules / cm^3]
    met.data["M"] = met.data["M"].transpose("latitude", "longitude", "level", "time")

    # Calculate H2O number concentration to feed into box model calcs
    met.data["H2O"] = (
        (met["specific_humidity"].data / constants.M_v) * N_A * rho_d * 1e-6
    )  # [molecules / cm^3]

    # Calculate O2 and N2 number concs based on M
    met.data["O2"] = 2.079e-01 * met.data["M"]
    met.data["N2"] = 7.809e-01 * met.data["M"]

    return met

# grab bg chem data
def grab_bg_chem(met):
    """Grab the background chemistry data for the month of the met data."""
    month = met["time"].data[0].dt.month

    bg_chem = xr.open_dataarray(
        os.environ['PYCONTRAILSDIR'] + "/models/boxmodel/species.nc"
    ).sel(month=month - 1)

    for s in [1, 2, 3, 5, 7, 9, 10, 13, 15, 16, 17, 18, 19, 20, 22, 24, 26, 27, 29, 31, 33, 35, 36, 37, 38, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 62, 63, 65, 66, 68, 69, 70, 72, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 200, 201, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219]:
        bg_chem[:, :, :, s-1] = 0


    bg_chem = bg_chem * 1e09  # convert mixing ratio to ppb

    return bg_chem

# convert latitude to latbox
def latitude_to_latbox(latitude):
        # Map the latitude to the range 0-1
        normalized_latitude = (latitude + 87.5) / 180

        # Map the normalized latitude to the range 1-72
        latbox = normalized_latitude * 36 + 1

        # Round to the nearest integer and return
        return round(latbox)

# convert longitude to longbox
def longitude_to_longbox(longitude):
        # Map the longitude to the range 0-1
        normalized_longitude = (longitude + 177.5) / 360

        # Map the normalized longitude to the range 1-144
        longbox = normalized_longitude * 72 + 1

        # Round to the nearest integer and return
        return round(longbox)

# convert altitude to pressure level
def get_pressure_level(alt):
        # Convert alt to pressure level (hPa)``
        chem_pressure_levels = np.array([962, 861, 759, 658, 556, 454, 353, 251, 150.5])

        # Convert altitude to pressure using a standard atmosphere model
        pressure = units.m_to_pl(alt)

        # Find the index of the closest value in the array
        idx = (np.abs(chem_pressure_levels - pressure)).argmin()

        return idx
