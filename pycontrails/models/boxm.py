"""A photochemical trajectory model for simulating atmospheric chemistry and NOx emissions."""

from __future__ import annotations

import os
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr
import dask as da
from boxm_f2py import boxm_f2py

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
from pycontrails.physics import constants, geo, units


@dataclass
class BoxmParams(ModelParams):
    """Default box model parameters."""

    lat_bounds: tuple[float, float] | None = None  # latmin, latmax [deg]
    lon_bounds: tuple[float, float] | None = None  # lonmin, lonmax [deg]
    alt_bounds: tuple[float, float] | None = None  # altmin, altmax [m]
    hres_chem: float = 1.0  # chemistry horizontal resolution [deg]
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
    bg_chem: MetDataset
    met_required = True

    def __init__(
        self,
        met: MetDataset,
        bg_chem: MetDataset,
        params: dict[str, Any] | None = None,
        **params_kwargs: Any,
    ) -> None:

        # call Model init
        super().__init__(met, params=params, **params_kwargs)

        # check if met and chem datasets are provided
        self.met = met
        self.bg_chem = bg_chem

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
        self.update_params(params)
        self.set_source(emi)
        self.source = self.require_source_type(MetDataset)

        self.process_datasets()
        self.to_csvs()
        self.run_boxm()

    def process_datasets(self):
        """Process the met, bg_chem, and emi datasets to prepare them for the box model."""
        # create lists for lats, lons, alts, and times
        lats_chem = np.arange(
            self.params.lat_bounds[0], self.params.lat_bounds[1], self.params.hres_chem
        )

        lons_chem = np.arange(
            self.params.lon_bounds[0], self.params.lon_bounds[1], self.params.hres_chem
        )

        alts_chem = np.arange(
            self.params.alt_bounds[0], self.params.alt_bounds[1], self.params.vres_chem
        )

        times_chem = pd.date_range(
            start=self.params.t0_chem,
            end=self.params.t0_chem + self.params.rt_chem,
            freq=self.params.ts_chem,
        )

        self.chem_dict = dict(lats=lats_chem, lons=lons_chem, alts=alts_chem, times=times_chem)

        # chunk met and emi data
        self.met.data = self.met.data.chunk(
            {"longitude": "auto", "latitude": "auto", "level": "auto", "time": "auto"}
        )

        self.source.data = self.source.data.chunk(
            {"longitude": "auto", "latitude": "auto", "level": "auto", "time": "auto"}
        )

        # interpolate and downselect met to chem grid
        self.met.data = self.met.data.interp(
            longitude=lons_chem,
            latitude=lats_chem,
            level=units.m_to_pl(alts_chem),
            time=times_chem,
            method="linear",
        )

        self.met = self.met.downselect(
            self.params.lon_bounds[0],
            self.params.lat_bounds[0],
            self.params.alt_bounds[0],
            self.params.lon_bounds[1],
            self.params.lat_bounds[1],
            self.params.alt_bounds[1],
        )

        # interpolate and downselect bg_chem to chem grid
        self.bg_chem.data = self.bg_chem.data.interp(
            longitude=lons_chem, latitude=lats_chem, level=units.m_to_pl(alts_chem), method="linear"
        )

        self.bg_chem = self.bg_chem.downselect(
            self.params.lon_bounds[0],
            self.params.lat_bounds[0],
            self.params.alt_bounds[0],
            self.params.lon_bounds[1],
            self.params.lat_bounds[1],
            self.params.alt_bounds[1],
        )

        # interpolate and downselect emi to chem grid
        self.source.data = self.source.data.interp(
            longitude=lons_chem,
            latitude=lats_chem,
            level=units.m_to_pl(alts_chem),
            time=times_chem,
            method="linear",
        )

        self.source = self.source.downselect(
            self.params.lon_bounds[0],
            self.params.lat_bounds[0],
            self.params.alt_bounds[0],
            self.params.lon_bounds[1],
            self.params.lat_bounds[1],
            self.params.alt_bounds[1],
        )

    def to_csvs(self):
        """Convert the met, bg_chem, and emi datasets to csv files for use in the box model."""
        # send dfs to csvs for chem analysis
        met_df = self.met.data.to_dask_dataframe(
            dim_order=["time", "level", "longitude", "latitude"]
        )
        met_df["latitude"] = met_df["latitude"].map("{:+08.3f}".format)
        met_df["longitude"] = met_df["longitude"].map("{:+08.3f}".format)
        met_df["sza"] = met_df["sza"].map("{:+0.3e}".format)
        met_df = met_df.apply(
            (lambda x: x.map("{:0.3e}".format) if x.dtype in ["float32", "float64"] else x), axis=1
        )

        bg_chem_df = self.bg_chem.data.to_dask_dataframe(
            dim_order=["level", "longitude", "latitude"]
        )
        bg_chem_df["latitude"] = bg_chem_df["latitude"].map("{:+08.3f}".format)
        bg_chem_df["longitude"] = bg_chem_df["longitude"].map("{:+08.3f}".format)
        bg_chem_df["month"] = bg_chem_df["month"].map("{:02d}".format)
        bg_chem_df = bg_chem_df.apply(
            (lambda x: x.map("{:0.3e}".format) if x.dtype in ["float32", "float64"] else x), axis=1
        )

        emi_df = self.source.data.to_dask_dataframe(
            dim_order=["time", "longitude", "latitude"]
        ).fillna(0)
        emi_df["latitude"] = emi_df["latitude"].map("{:+08.3f}".format)
        emi_df["longitude"] = emi_df["longitude"].map("{:+08.3f}".format)
        emi_df = emi_df.apply(
            (lambda x: x.map("{:0.3e}".format) if x.dtype in ["float32", "float64"] else x), axis=1
        )

        # Remove temporary files if they exist
        if os.path.exists("met_df.csv"):
            os.remove("met_df.csv")
        if os.path.exists("bg_chem_df.csv"):
            os.remove("bg_chem_df.csv")
        if os.path.exists("emi_df.csv"):
            os.remove("emi_df.csv")

        # Write DataFrame 1 to the temporary file
        met_df.to_csv("met_df.csv", single_file=True, index=False)

        # Write DataFrame 2 to the temporary file
        bg_chem_df.to_csv("bg_chem_df.csv", single_file=True, index=False)

        # Write DataFrame 3 to the temporary file
        emi_df.to_csv("emi_df.csv", single_file=True, index=False)

    def run_boxm(self):
        """Run the box model in fortran using f2py interface."""
        ncell = (
            len(self.chem_dict["lats"]) * len(self.chem_dict["lons"]) * len(self.chem_dict["alts"])
        )

        nts = len(self.chem_dict["times"]) - 1

        dts = self.params.ts_chem.total_seconds()

        boxm_f2py.init(ncell)

        for t in range(nts):
            # if t*dts % 3600 == 0:
            print("Time: ", t)

            boxm_f2py.read(ncell)
            boxm_f2py.calc_aerosol()
            boxm_f2py.chemco()
            boxm_f2py.calc_j(ncell)
            boxm_f2py.photol()

            if t != 0:
                boxm_f2py.deriv(dts)

            boxm_f2py.write(dts, ncell)

        boxm_f2py.deallocate()


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
    return