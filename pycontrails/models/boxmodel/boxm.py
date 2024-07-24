"""A photochemical trajectory model for simulating atmospheric chemistry and NOx emissions."""

from __future__ import annotations

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


@dataclass
class BoxmParams(ModelParams):
    """Default box model parameters."""

    lat_bounds: tuple[float, float] | None = None  # latmin, latmax [deg]
    lon_bounds: tuple[float, float] | None = None  # lonmin, lonmax [deg]
    alt_bounds: tuple[float, float] | None = None  # altmin, altmax [m]
    hres_chem: float = 0.5  # chemistry horizontal resolution [deg]
    vres_chem: float = 500  # chemistry vertical resolution [m]

    t0_chem: datetime | None = None  # chemistry start time
    rt_chem: timedelta | None = None  # chemistry runtime
    ts_chem: timedelta | None = None  # chemistry timestep


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
        # bg_chem: MetDataset,
        params: dict[str, Any] | None = None,
        **params_kwargs: Any,
    ) -> None:

        # call Model init
        super().__init__(met, params=params, **params_kwargs)

        # check if met and chem datasets are provided
        if met:
            self.met = met

        # if bg_chem:
        #     self.bg_chem = bg_chem
        # else:
        self.bg_chem = grab_bg_chem(met)

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
        self.update_params(params)
        self.set_source(emi)
        self.source = self.require_source_type(MetDataset)

        self.process_datasets()
        self.init_boxm_ds()
        self.to_netcdfs()
        self.run_boxm()
        self.run_boxm_orig()


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

        self.boxm_ds["Y"] = (["time", "level", "longitude", "latitude", "species"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], self.boxm_ds.dims["species"])))

        self.boxm_ds["Y_orig"] = (["time", "level", "longitude", "latitude", "species"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], self.boxm_ds.dims["species"])))

        self.boxm_ds["FL"] = (["time", "level", "longitude", "latitude", "flux_rates"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], 130)))

    def to_netcdfs(self):
        """Convert the met, bg_chem, and emi datasets to boxm_ds.nc for use in the box model."""

        # stack datasets to get cell index for fortran
        self.boxm_ds_stacked = self.boxm_ds.stack(
            {"cell": ["level", "longitude", "latitude"]}
        )
        self.boxm_ds_stacked = self.boxm_ds_stacked.reset_index("cell")

        print(self.boxm_ds_stacked)

        # Delete any existing netCDF files
        if pathlib.Path(self.path + "boxm_ds.nc").exists():
            print("deleting boxm_ds.nc")
            pathlib.Path(self.path + "boxm_ds.nc").unlink()

        if pathlib.Path(self.path + "boxm_input.txt").exists():
            print("deleting boxm_input.txt")
            pathlib.Path(self.path + "boxm_input.txt").unlink()

            # delete any input and output files
        if pathlib.Path(self.path + "BACKITNE.OUT").exists():
            print("deleting BACKITNE.OUT")
            pathlib.Path(self.path + "BACKITNE.OUT").unlink()

        if pathlib.Path(self.path + "ZEN.OUT").exists():
            print("deleting ZEN.OUT")
            pathlib.Path(self.path + "ZEN.OUT").unlink()

        if pathlib.Path(self.path + "J.OUT").exists():
            print("deleting J.OUT")
            pathlib.Path(self.path + "J.OUT").unlink()

        if pathlib.Path(self.path + "DJ.OUT").exists():
            print("deleting DJ.OUT")
            pathlib.Path(self.path + "DJ.OUT").unlink()

        if pathlib.Path(self.path + "RC.OUT").exists():
            print("deleting RC.OUT")
            pathlib.Path(self.path + "RC.OUT").unlink()

        if pathlib.Path(self.path + "Y.OUT").exists():
            print("deleting Y.OUT")
            pathlib.Path(self.path + "Y.OUT").unlink()

        if pathlib.Path(self.path + "FL.OUT").exists():
            print("deleting FL.OUT")
            pathlib.Path(self.path + "FL.OUT").unlink()

        # Convert DataFrames to Datasets and write to netCDF
        self.boxm_ds_stacked.to_netcdf(self.path + "boxm_ds.nc", mode="w")

    def run_boxm(self):
        """Run the box model in fortran using subprocess."""
        
        subprocess.call(
            [self.path + "boxm"]
        )
        
        # open nc file
        self.boxm_ds = xr.open_dataset(self.path + "boxm_ds.nc")
        print(self.boxm_ds)

    def run_boxm_orig(self):
        """Run the original box model in fortran using subprocess."""
        
        for cell in range(1):
        


            # generate input file
            self.gen_boxm_input_file(cell)

            # calls fortran with input file, and generates .OUT files
            subprocess.call(
                [self.path + "boxm_orig"]
            )

            # read output files to nc
            self.output_to_nc(cell)

        # convert lat lon to coords
        # chem = chem.set_coords(("level", "longitude", "latitude"))

        # chem = chem.assign_coords(
        #     time=self.met["time"].data.values[: chem.dims["time"]],
        #     level=chem.level,
        #     longitude=chem.longitude,
        #     latitude=chem.latitude,
        #     cell=range(),
            # "level": self.met["level"].data.values,
            # "longitude": self.met["longitude"].data.values,
            # "latitude": self.met["latitude"].data.values}
  

        # chem = chem.set_index(cell=["level", "longitude", "latitude"])
        # chem = chem.unstack("cell")

        # chem = chem.set_coords(['level', 'longitude', 'latitude'])


    def gen_boxm_input_file(self, cell):
        """Generate the input file for the original box model."""

        # delete any existing input files
        if pathlib.Path(self.path + "boxm_input.txt").exists():
                pathlib.Path("/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/boxm_input.txt").unlink()

        # open file
        boxm_input = open(self.path + "boxm_input.txt", "w")

        start_time = pd.to_datetime(self.met["time"].data.values[0])
        end_time = pd.to_datetime(self.met["time"].data.values[-1])
        runtime = int((end_time - start_time) / np.timedelta64(1, 'D')) % 365
        day = start_time.day
        month = start_time.month
        year = start_time.year
        altitude = self.boxm_ds_stacked["altitude"].values[cell]
        plevel = self.boxm_ds_stacked["level"].values[cell]
        level = get_pressure_level(altitude)
        longitude = self.boxm_ds_stacked["longitude"].values[cell]
        longbox = longitude_to_longbox(longitude)
        latitude = self.boxm_ds_stacked["latitude"].values[cell]
        latbox = latitude_to_latbox(latitude)
        M = self.boxm_ds_stacked["M"].values[0,cell]
        P = self.boxm_ds_stacked["air_pressure"].values[cell]
        H2O = self.boxm_ds_stacked["H2O"].values[0,cell]
        temp = self.boxm_ds_stacked["air_temperature"].values[0,cell]

        boxm_input.write(repr(day) + "\n" + repr(month) + "\n" + repr(year) + "\n" + repr(level) + "\n" + repr(longbox) + "\n" + repr(latbox) + "\n" + repr(runtime) + "\n" + repr(M) + "\n" + repr(plevel) + "\n" + repr(H2O) + "\n" + repr(temp) + "\n")

        for s in ["NO2", "NO", "O3", "CO", "CH4", "HCHO", "CH3CHO", "CH3COCH3",
                           "C2H6", "C2H4", "C3H8", "C3H6", "C2H2", "NC4H10", "TBUT2ENE",
                           "BENZENE", "TOLUENE", "OXYL", "C5H8", "H2O2", "HNO3", "C2H5CHO",
                           "CH3OH", "MEK", "CH3OOH", "PAN", "MPAN", "OH"]:
            
            boxm_input.write(repr(self.bg_chem.sel(species=s).loc[latitude, 
                                                           longitude, plevel].values.item()) + "\n")
            
        boxm_input.close()

    
    def output_to_nc(self, cell):
         
        # ZEN_df = pd.read_csv("/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/ZEN.OUT", header=0,
        #                   names=['TIME', 'ZEN', 'TEMP', 'cosx', 'secx'], dtype=np.float64)
         
        # J_df = pd.read_csv("/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/J.OUT", header=0,
        #                     names=['TIME', 'J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'J11', 'J12', 'J13','J14', 'J15', 'J16', 'J17', 'J18', 'J19', 'J20', 'J21', 'J22', 'J23', 'J24', 'J25', 'J26', 'J27', 'J28', 'J29', 'J30', 'J31', 'J32', 'J33', 'J34', 'J35', 'J36', 'J37', 'J38', 'J39', 'J40', 'J41', 'J42', 'J43', 'J44', 'J45', 'J46', 'J47', 'J48', 'J49', 'J50'], dtype=np.float64)

        # DJ_df = pd.read_csv("/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/DJ.OUT", header=0,
        #                         names=['TIME', 'DJ1', 'DJ2', 'DJ3', 'DJ4', 'DJ5', 'DJ6', 'DJ7', 'DJ8', 'DJ9', 'DJ10', 'DJ11', 'DJ12', 'DJ13','DJ14', 'DJ15', 'DJ16', 'DJ17', 'DJ18', 'DJ19', 'DJ20', 'DJ21', 'DJ22', 'DJ23', 'DJ24', 'DJ25', 'DJ26', 'DJ27', 'DJ28', 'DJ29', 'DJ30', 'DJ31', 'DJ32', 'DJ33', 'DJ34', 'DJ35', 'DJ36', 'DJ37', 'DJ38', 'DJ39', 'DJ40', 'DJ41', 'DJ42', 'DJ43', 'DJ44', 'DJ45', 'DJ46', 'DJ47', 'DJ48', 'DJ49', 'DJ50'], dtype=np.float64)

        # RC_df = pd.read_csv("/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/RC.OUT", header=0,
        #                     names=['TIME', 'RC1', 'RC2', 'RC3', 'RC4', 'RC5', 'RC6', 'RC7', 'RC8', 'RC9', 'RC10', 'RC11', 'RC12', 'RC13','RC14', 'RC15', 'RC16', 'RC17', 'RC18', 'RC19', 'RC20', 'RC21', 'RC22', 'RC23', 'RC24', 'RC25', 'RC26', 'RC27', 'RC28', 'RC29', 'RC30', 'RC31', 'RC32', 'RC33', 'RC34', 'RC35', 'RC36', 'RC37', 'RC38', 'RC39', 'RC40', 'RC41', 'RC42', 'RC43', 'RC44', 'RC45', 'RC46', 'RC47', 'RC48', 'RC49', 'RC50'], dtype=np.float64)

        Y_df = pd.read_csv("/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/Y.OUT", header=0,
                            names=['TIME', 'O1D', 'O', 'OH', 'NO2', 'NO3', 'O3', 'N2O5', 'NO', 'HO2', 'H2', 'CO', 'H2O2', 'HONO', 'HNO3', 'HO2NO2', 'SO2', 'SO3', 'HSO3', 'NA', 'SA', 'CH4', 'CH3O2', 'C2H6', 'C2H5O2', 'C3H8', 'IC3H7O2', 'RN10O2', 'NC4H10', 'RN13O2', 'C2H4', 'HOCH2CH2O2', 'C3H6', 'RN902', 'TBUT2ENE', 'RN12O2', 'NRN6O2', 'NRN9O2', 'NRN12O2', 'HCHO', 'HCOOH', 'CH3CO2H', 'CH3CHO', 'C5H8', 'RU14O2', 'NRU14O2',
                            'UCARB10', 'APINENE', 'RTN28O2', 'NRTN28O2', 'RTN26O2'], dtype=np.float64)
        print(Y_df["TIME"].values)
        

        for s, species in enumerate(Y_df.columns[1:]):

            self.boxm_ds["Y_orig"].loc[:, species, cell] = Y_df[species].values



    # animate chemdataset
    def anim_chem(self, mda):
        """Animate the chemical concentrations."""
        fig, (ax, cbar_ax) = plt.subplots(
            1, 2, gridspec_kw={"width_ratios": (0.9, 0.05), "wspace": 0.2}, figsize=(12, 8)
        )

        times = mda["time"].values

        def heatmap_func(t):
            ax.cla()
            ax.set_title(t)

            mda.sel(time=t).transpose("latitude", "longitude").plot(
                ax=ax, cbar_ax=cbar_ax, add_colorbar=True, vmin=mda.min(), vmax=mda.max()
            )

        anim = FuncAnimation(fig, heatmap_func, frames=times, blit=False)

        filename = pathlib.Path("plume.gif")

        anim.save(filename, dpi=300, writer=PillowWriter(fps=8))


# functions used in boxm
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
        "/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/species.nc"
    ).sel(month=month - 1)

    bg_chem = bg_chem * 1e09  # convert mixing ratio to ppb

    return bg_chem


def latitude_to_latbox(latitude):
        # Map the latitude to the range 0-1
        normalized_latitude = (latitude + 87.5) / 180

        # Map the normalized latitude to the range 1-72
        latbox = normalized_latitude * 36

        # Round to the nearest integer and return
        return round(latbox)

def longitude_to_longbox(longitude):
        # Map the longitude to the range 0-1
        normalized_longitude = (longitude + 177.5) / 360

        # Map the normalized longitude to the range 1-144
        longbox = normalized_longitude * 72

        # Round to the nearest integer and return
        return round(longbox)

def get_pressure_level(alt):
        # Convert alt to pressure level (hPa)``
        chem_pressure_levels = np.array([962, 861, 759, 658, 556, 454, 353, 251, 150.5])

        # Convert altitude to pressure using a standard atmosphere model
        pressure = units.m_to_pl(alt)

        # Find the index of the closest value in the array
        idx = (np.abs(chem_pressure_levels - pressure)).argmin()

        return idx