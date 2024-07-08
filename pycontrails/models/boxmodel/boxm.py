"""A photochemical trajectory model for simulating atmospheric chemistry and NOx emissions."""

from __future__ import annotations

import os
import pathlib
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.animation import FuncAnimation, PillowWriter

from pycontrails.core.met import MetDataArray, MetDataset
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
        self.update_params(params)
        self.set_source(emi)
        self.source = self.require_source_type(MetDataset)

        self.process_datasets()
        self.to_netcdfs()
        #chem = self.run_boxm()

        #return chem

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
            calc_sza(met["latitude"].data.values, met["longitude"].data.values, met["time"].data.values),
        )

        # interpolate and downselect bg_chem to chem grid
        self.bg_chem = self.bg_chem.interp(
            longitude=met["longitude"].data,
            latitude=met["latitude"].data,
            level=met["level"].data,
            method="linear",
        )

        # # interpolate and downselect emi to chem grid
        # self.source.data = self.source.data.interp(
        #     longitude=met["longitude"].data,
        #     latitude=met["latitude"].data,
        #     time=met["time"].data,
        #     method="linear",
        # )

    def to_netcdfs(self):
        """Convert the met, bg_chem, and emi datasets to boxm_input.nc for use in the box model."""

        path = "/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/"
        self.boxm_input = xr.merge([self.met.data, self.bg_chem, self.source.data])
        self.boxm_input = self.boxm_input.drop_vars(["specific_humidity", "relative_humidity", "eastward_wind", "northward_wind", "lagrangian_tendency_of_air_pressure", "month"])

        # stack datasets to get cell index for fortran
        self.boxm_input_stacked = self.boxm_input.stack({"cell": ["level", "longitude", "latitude"]})
        self.boxm_input_stacked = self.boxm_input_stacked.reset_index("cell")

        # Convert DataFrames to Datasets and write to netCDF
        self.boxm_input_stacked.to_netcdf(path + "boxm_input.nc", mode="w")

    def run_boxm(self):
        """Run the box model in fortran using f2py interface."""
        ncell = len(self.met["latitude"]) * len(self.met["longitude"]) * len(self.met["altitude"])

        nts = len(self.met["time"]) - 1

        dts = self.params["ts_chem"].total_seconds()

        print(ncell, nts, dts)

        boxm.boxm.init(ncell, nts)

        for ts in range(nts):
            # if t*dts % 3600 == 0:
            print("Time: ", self.met["time"].data[ts].values)

            boxm.boxm.read()
            boxm.boxm.calc_aerosol()
            boxm.boxm.chemco()
            boxm.boxm.calc_j()
            boxm.boxm.photol()

            if ts != 0:
                boxm.boxm.deriv()

            boxm.boxm.write()

        boxm.boxm.deallocate()

        # open chem dataset
        chem = xr.open_dataset("/home/ktait98/pycontrails_kt/pycontrails/models/boxmodel/chem.nc")

        # convert lat lon to coords
        # chem = chem.set_coords(("level", "longitude", "latitude"))

        chem = chem.assign_coords(
            time=self.met["time"].data.values[: chem.dims["time"]],
            level=chem.level,
            longitude=chem.longitude,
            latitude=chem.latitude,
            cell=range(ncell)
            # "level": self.met["level"].data.values,
            # "longitude": self.met["longitude"].data.values,
            # "latitude": self.met["latitude"].data.values}
        )

        chem = chem.set_index(cell=["level", "longitude", "latitude"])
        chem = chem.unstack("cell")

        # chem = chem.set_coords(['level', 'longitude', 'latitude'])

        return chem


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


    
