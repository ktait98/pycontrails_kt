"""Gridded Plume Analysis Tool (GPAT).

Simulate aircraft trajectories, estimate aircraft performance, fuel burn and emissions.

Plot associated aircraft exhaust plumes, subject to Gaussian dispersion and advection.
Aggregate plumes to an Eulerian grid for photochemical and microphysical processing.
"""

import argparse
import os
import pathlib
import pickle
import random
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field, fields
from typing import Literal, Optional

import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Geod
import json

from pycontrails.core import Flight, GeoVectorDataset, MetDataset, models
from pycontrails.core.models import Model
from pycontrails.models.dry_advection import DryAdvection
from pycontrails.models.emissions import Emissions
from pycontrails.models.gpat.plume_to_grid import plume_to_grid
from pycontrails.models.ps_model import PSFlight
from pycontrails.physics import constants, geo, thermo, units


### GPAT Model Parameters ###
@dataclass
class SimParams:
    """Default simulation parameters."""

    # Temporal domain
    # flight time
    t_fl: tuple[pd.Timestamp, pd.Timedelta, pd.Timedelta] = field(
        default_factory=lambda: (
            pd.to_datetime("2022-01-20 13:00:00"),
            pd.Timedelta(minutes=2),
            pd.Timedelta(hours=1),
        )
    )  # (start time, time step, run time)

    # plume time
    t_pl: tuple[pd.Timestamp, pd.Timedelta, pd.Timedelta] = field(
        default_factory=lambda: (
            pd.to_datetime("2022-01-20 13:00:00"),
            pd.Timedelta(minutes=2),
            pd.Timedelta(hours=2),
        )
    )  # (start time, time step, max age)

    # simulation time
    t_sim: tuple[pd.Timestamp, pd.Timedelta, pd.Timedelta] = field(
        default_factory=lambda: (
            pd.to_datetime("2022-01-20 12:00:00"),
            pd.Timedelta(seconds=20),
            pd.Timedelta(hours=120),
        )
    )  # (start time, time step, run time)

    #  spatial domain
    lat_bounds: tuple[float, float] = (0.0, 1.0)  # lat bounds [deg]
    lon_bounds: tuple[float, float] = (0.0, 1.0)  # lon bounds [deg]
    alt_bounds: tuple[float, float] = (12000, 13000)  # alt bounds [m]
    hres_sim_c: float = 0.01  # horizontal resolution [deg]
    vres_sim_c: float = 500  # vertical resolution [m]
    hres_sim_f: float = 0.001  # horizontal resolution [deg]
    vres_sim_f: float = 100  # vertical resolution [m]

    run_path: Optional[str] = None  # path to run GPAT from
    data_path: Optional[str] = None  # path to data directory
    job_id: Optional[str] = None  # job ID
    
    # run_chem: bool = True  # whether to run chemistry model
    # wind_effects: bool = True  # whether to include wind effects in plume dispersion

@dataclass
class FlParams:
    """Default flight/fleet parameters."""

    mode: Literal["direct", "synthetic"] = "direct"
    file: Optional[str] = None

    ac_type: Optional[str] = "A320"  # aircraft type
    fl0_speed: Optional[float] = 100.0  # m/s
    fl0_heading: Optional[float] = 0.0  # deg
    fl0_coords0: Optional[tuple[float, float, float]] = (0.1, 0.125, 12500)  # lat, lon, alt [deg, deg, m]
    sep_dist: Optional[tuple[float, float, float]] = (5000, 2000, 0)  # dx, dy, dz [m]
    n_ac: Optional[int] = 1  # number of aircraft

@dataclass
class PlumeParams:
    """Default plume dispersion parameters."""

    depth: float = 50.0  # initial plume depth, [m]
    width: float = 50.0  # initial plume width, [m]
    verbose_outputs: bool = False  # print verbose outputs
    n_slices: int = 5  # number of slices in the plume
    shear: float = 0.01  # shear [m/s]

@dataclass
class MetParams:
    """Default meteorological parameters."""

    eastward_wind: Optional[float] = 0.0  # m/s
    northward_wind: Optional[float] = 0.0  # m/s
    lagrangian_tendency_of_air_pressure: Optional[float] = 0.0  # Pa/s

@dataclass
class ChemParams:
    """Default chemistry parameters."""
    run_chem: bool = True  # whether to run chemistry model
    species_in: tuple = ("NO",)
    species_out: tuple = ("O3", "NO2", "NO", "NO3", "N2O5", "HNO3", 
                          "HONO", "HO2", "OH", "H2O2", "CO", "CH4", "CH3O2")

# @dataclass
# class ContrailParams:

class GPAT(Model):
    """Gridded Plume Analysis Tool (GPAT).

    Simulate aircraft trajectories, estimate aircraft performance, fuel burn and emissions. Then
    aggregates emissions, bg chemistry and meteorology to an Eulerian grid for photochemical and
    microphysical processing.

    Parameters
    ----------
    sim_params : SimParams
        Simulation parameters.
    fl_params : FlParams
        Flight parameters.
    plume_params : PlumeParams
        Plume dispersion parameters.
    met_params : MetParams
        Meteorological parameters.
    chem_params : ChemParams
        Chemistry parameters.
    contrail_params : ContrailParams
        Contrail parameters.
    """

    name = "GPAT"
    long_name = "Gridded Plume Analysis Tool"
    # default_params = (FlParams, PlumeParams, SimParams)

    def __init__(self, sim_params: SimParams, 
                 fl_params: FlParams,
                 plume_params: PlumeParams, 
                 met_params: MetParams, 
                 chem_params: ChemParams):
        super().__init__()

        # Generate the coarse grid vectors
        self.lats = np.arange(
            sim_params.lat_bounds[0] + sim_params.hres_sim_c / 2,
            sim_params.lat_bounds[1],
            sim_params.hres_sim_c,
        )
        self.lons = np.arange(
            sim_params.lon_bounds[0] + sim_params.hres_sim_c / 2,
            sim_params.lon_bounds[1],
            sim_params.hres_sim_c,
        )
        self.alts = np.arange(
            sim_params.alt_bounds[0] + sim_params.vres_sim_c / 2,
            sim_params.alt_bounds[1],
            sim_params.vres_sim_c,
        )

        self.levels = units.m_to_pl(self.alts)

        # Generate time vectors
        self.times_fl = pd.date_range(
            start=sim_params.t_fl[0],
            end=sim_params.t_fl[0] + sim_params.t_fl[2],
            freq=sim_params.t_fl[1],
        )

        self.times_pl = pd.date_range(
            start=sim_params.t_pl[0],
            end=sim_params.t_pl[0] + sim_params.t_pl[2],
            freq=sim_params.t_pl[1],
        )

        self.times_sim = pd.date_range(
            start=sim_params.t_sim[0],
            end=sim_params.t_sim[0] + sim_params.t_sim[2],
            freq=sim_params.t_sim[1],
        )
        
        if sim_params.run_path is None:
            self.run_path = os.environ["PYCONTRAILSDIR"] + "models/gpat/"

        else:
            self.run_path = sim_params.run_path

        if sim_params.data_path is None:
            self.data_path = "/projects/Impact_of_aviation_on_climate/Kieran2024/"

        else:
            self.data_path = sim_params.data_path

        if sim_params.job_id is None:
            try:
                self.job_id = os.environ["SLURM_JOB_ID"]
            except KeyError:
                # If SLURM_JOB_ID is not found, generate a random number as job ID
                self.job_id = str(random.randint(100000, 999999))

        else:
            self.job_id = sim_params.job_id

        sim_params.date_created = pd.Timestamp.now()
        chem_params.species_out_num = grab_species_num(self.run_path, chem_params.species_out)
        
        self.inputs_job = self.data_path + "inputs/" + self.job_id + "/"
        self.inputs_glob = self.data_path + "inputs/glob/"
        self.outputs_job = self.data_path + "outputs/" + self.job_id + "/"  

        if os.path.exists(self.inputs_job):
            shutil.rmtree(self.inputs_job)

        if os.path.exists(self.outputs_job):
            shutil.rmtree(self.outputs_job)

        os.makedirs(self.inputs_job)
        os.makedirs(self.outputs_job)      

        all_params = {
            "sim_params": sim_params,
            "fl_params": fl_params,
            "plume_params": plume_params,
            "met_params": met_params,
            "chem_params": chem_params,
        }

        # Set the model parameters
        self.sim_params = sim_params
        self.fl_params = fl_params
        self.plume_params = plume_params
        self.met_params = met_params
        self.chem_params = chem_params
        self.all_params = all_params

    def eval(self):
        """Run the GPAT model."""

        # Generate flight trajectory points
        self.fl = self.traj_gen()

        # Generate meteorological data
        self.met = self.gen_met()

        # Generate background chemistry data
        self.bg_chem = self.gen_bg_chem()

        # Calculate aircraft performance using PS Model
        self.fl = self.ac_perf()

        # Estimate emissions using Pycontrails Emissions Model
        self.fl = self.emissions()

        # Simulate plume dispersion/advection using Pycontrails Dry Advection Model
        self.fl, self.pl = self.sim_plumes()

        # Run BOXM
        self.chem = self.run_boxm()

        # Generate outputs
        # self.gen_outputs()

    # Model methods
    def traj_gen(self) -> list[Flight]:
        """Generate flight trajectory points."""
        fl_params = self.fl_params
        sim_params = self.sim_params

        if fl_params.mode == "direct":
            # provide flight file - csv to convert to Flight object (Pd df)
            if fl_params.file is None:
                raise ValueError("Flight file must be provided for direct mode.")
            fl = Flight.read_csv(fl_params.file)
            fl.attrs = {"flight_id": 0, "aircraft_type": fl_params.ac_type}
            mask = (
                (fl["latitude"] > self.sim_params.lat_bounds[0] + 0.01)
                & (fl["latitude"] < self.sim_params.lat_bounds[1] - 0.01)
                & (fl["longitude"] > self.sim_params.lon_bounds[0] + 0.01)
                & (fl["longitude"] < self.sim_params.lon_bounds[1] - 0.01)
                & (fl["altitude"] > self.sim_params.alt_bounds[0])
                & (fl["altitude"] < self.sim_params.alt_bounds[1])
            )
            fl = fl.filter(mask)
            fl = [fl]
            return fl

        # generate synthetic formation flight
        if fl_params.mode == "synthetic":
            if fl_params.n_ac < 1:
                raise ValueError("Number of aircraft must be at least 1.")
            if fl_params.ac_type is None:
                raise ValueError("Aircraft type must be provided for synthetic mode.")
            if fl_params.fl0_coords0 is None:
                raise ValueError("Initial coordinates must be provided for synthetic mode.")
            if fl_params.fl0_speed is None:
                raise ValueError("Flight speed must be provided for synthetic mode.")
            if fl_params.fl0_heading is None:
                raise ValueError("Flight heading must be provided for synthetic mode.")
            if fl_params.sep_dist is None:
                raise ValueError("Separation distances must be provided for synthetic mode.")
            
            fl = []

            lat0, lon0, alt0 = fl_params.fl0_coords0
            heading = fl_params.fl0_heading
            dist = fl_params.fl0_speed * sim_params.t_fl[2].total_seconds()

            # calculate the final coordinates
            geod = Geod(ellps="WGS84")
            lon1, lat1, _ = geod.fwd(lon0, lat0, heading, dist)

            # create flight object for leader flight and resample points according to ts_fl
            df = pd.DataFrame()
            df["longitude"] = [lon0, lon1]
            df["latitude"] = [lat0, lat1]
            df["altitude"] = [alt0, alt0]
            df["time"] = [sim_params.t_fl[0], (sim_params.t_fl[0] + sim_params.t_fl[2])]

            ts_fl_min = int(sim_params.t_fl[1].total_seconds() / 60)

            fl0 = Flight(df).resample_and_fill(freq=f"{ts_fl_min}min")
            fl0.attrs = {"flight_id": 0, "aircraft_type": fl_params.ac_type}
            mask = (
                (fl0["latitude"] > self.sim_params.lat_bounds[0] + 0.01)
                & (fl0["latitude"] < self.sim_params.lat_bounds[1] - 0.01)
                & (fl0["longitude"] > self.sim_params.lon_bounds[0] + 0.01)
                & (fl0["longitude"] < self.sim_params.lon_bounds[1] - 0.01)
                & (fl0["altitude"] > self.sim_params.alt_bounds[0])
                & (fl0["altitude"] < self.sim_params.alt_bounds[1])
            )

            fl0 = fl0.filter(mask)
            fl.append(fl0)

            fli = fl0

            if fl_params.n_ac > 1:
                # create follower flight trajectories
                for i in range(1, fl_params.n_ac):
                    fli = fli.copy()

                    # calculate new coords for follower flight
                    lon_dx, lat_dx, _ = geod.fwd(lon0, lat0, heading, fl_params.sep_dist[0])
                    lon_dx_dy, lat_dx_dy, _ = geod.fwd(
                        lon_dx, lat_dx, heading + 90, fl_params.sep_dist[1]
                    )
                    alt_dx_dy = alt0 + fl_params.sep_dist[2]

                    # Calculate the differences in lat, lon, alt
                    dlat = lat_dx_dy - lat0
                    dlon = lon_dx_dy - lon0
                    dalt = alt_dx_dy - alt0

                    # Update the latitude and longitude of each point in the flight path
                    fli["latitude"] += dlat
                    fli["longitude"] += dlon
                    fli["altitude"] += dalt
                    fli.attrs = {"flight_id": int(i), "aircraft_type": fl_params.ac_type}

                    mask = (
                        (fli["latitude"] > self.sim_params.lat_bounds[0] + 0.01)
                        & (fli["latitude"] < self.sim_params.lat_bounds[1] - 0.01)
                        & (fli["longitude"] > self.sim_params.lon_bounds[0] + 0.01)
                        & (fli["longitude"] < self.sim_params.lon_bounds[1] - 0.01)
                        & (fli["altitude"] > self.sim_params.alt_bounds[0])
                        & (fli["altitude"] < self.sim_params.alt_bounds[1])
                    )
                    fli = fli.filter(mask)
                    fl.append(fli)

                    # Update starting coordinates for next flight
                    lon0, lat0, alt0 = lon_dx_dy, lat_dx_dy, alt_dx_dy

            return fl

    def gen_met(self) -> MetDataset:
        """Generate meteorology data."""
        met_params = self.met_params

        # Step 1: Create with STANDARD names for MetDataset validation
        met_standard = xr.Dataset(
            data_vars={
                "eastward_wind": (
                    ("time", "level", "latitude", "longitude"),
                    np.full((len(self.times_sim), len(self.levels), len(self.lats), len(self.lons)), met_params.eastward_wind),
                ),
                "northward_wind": (
                    ("time", "level", "latitude", "longitude"),
                    np.full((len(self.times_sim), len(self.levels), len(self.lats), len(self.lons)), met_params.northward_wind),
                ),
                "lagrangian_tendency_of_air_pressure": (
                    ("time", "level", "latitude", "longitude"),
                    np.full((len(self.times_sim), len(self.levels), len(self.lats), len(self.lons)), met_params.lagrangian_tendency_of_air_pressure),
                ),
            },
            coords={
                "longitude": self.lons,
                "latitude": self.lats,
                "level": self.levels,
                "time": self.times_sim,
            },
        )

        # Step 2: Initialize MetDataset (validates standard names)
        met = MetDataset(met_standard)

        month = self.times_sim[0].month

        # Step 3: Load and interpolate climatology with standard names
        air_temperature = (
            xr.open_dataarray(self.inputs_glob + "air_temperature.nc", engine="netcdf4")
            .sel(month=month - 1)
            .interp(
                longitude=self.lons,
                latitude=self.lats,
                level=self.levels,
                method="linear"
            )
            .broadcast_like(met.data["eastward_wind"])
        )

        h2o_concs = (
            xr.open_dataarray(self.inputs_glob + "h2o_concs.nc", engine="netcdf4")
            .sel(month=month - 1)
            .interp(
                longitude=self.lons,
                latitude=self.lats,
                level=self.levels,
                method="linear"
            )
            .broadcast_like(met.data["eastward_wind"])
        )
        N_A = 6.022e23  # Avogadro's number
        
        # Add temp and H2O to met dataset
        met.data["air_temperature"] = air_temperature
        met.data["H2O"] = h2o_concs.transpose("latitude", "longitude", "level", "time")

        # Calculate specific humidity and relative humidity
        rho_d = met["air_pressure"].data / (constants.R_d * met["air_temperature"].data)
        met.data["specific_humidity"] = met.data["H2O"] * constants.M_d / (N_A * rho_d * 1e-6)
        met.data["relative_humidity"] = thermo.rhi(
            met.data["specific_humidity"], met.data["air_temperature"], met.data["air_pressure"]
        )

        # Calculate number density of air (M) to feed into box model calcs
        met.data["M"] = (N_A / constants.M_d) * rho_d * 1e-6  # [molecules / cm^3]
        met.data["M"] = met.data["M"].transpose("latitude", "longitude", "level", "time")

        # Calculate O2 and N2 number concs based on M
        met.data["O2"] = 2.079e-01 * met.data["M"]
        met.data["N2"] = 7.809e-01 * met.data["M"]

        # calculate solar zenith angle
        met.data["sza"] = (
            ("latitude", "longitude", "time"),
            calc_sza(
                met["latitude"].data.values, met["longitude"].data.values, met["time"].data.values
            ),
        )

        return met

    def gen_bg_chem(self) -> xr.Dataset:
        """Generate background chemistry data."""
        month = self.times_sim[0].month

        bg_chem = (
            xr.open_dataset(self.inputs_glob + "species.nc", engine="netcdf4")
            .sel(month=month - 1)
        )

        for s in [1,2,3,5,7,9,10,13,15,16,17,18,19,20,22,24,26,27,29,31,33,35,36,37,38,40,
                  41,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,60,62,63,65,66,68,69,70,
                  72,74,75,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,
                  98,99,100,102,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,
                  119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,
                  137,138,139,140,141,142,143,145,146,147,148,149,150,151,152,153,154,155,
                  156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,
                  174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
                  192,193,194,195,196,197,199,200,201,203,204,205,206,207,208,209,210,211,
                  212,213,214,215,216,217,218,219]:
            bg_chem.bg_chem[:, :, :, s - 1] = 0

        bg_chem = bg_chem * 1e09  # convert mixing ratio to ppb

        # downselect and interpolate bg_chem to the simulation grid
        return bg_chem.interp(longitude=self.lons, latitude=self.lats, level=self.levels)

    def ac_perf(self) -> list[Flight]:
        """Calculate aircraft performance using PS Model."""
        met = self.met
        fl = self.fl

        ps_model = PSFlight()

        for i, fli in enumerate(fl):
            # downselect met data to the flight trajectory
            fli.downselect_met(met)
            fl[i]["air_temperature"] = models.interpolate_met(met, fli, "air_temperature")
            fl[i]["specific_humidity"] = models.interpolate_met(met, fli, "specific_humidity")
            fl[i]["true_airspeed"] = fli.segment_groundspeed()
            print(f"flight {i} done")

            # get ac performance data using Poll-Schumann Model
            fl[i] = ps_model.eval(fl[i])

        return fl

    def emissions(self) -> list[Flight]:
        """Estimate emissions using Pycontrails Emissions Model."""
        sim_params = self.sim_params
        # met = self.met
        fl = self.fl

        emi_model = Emissions()

        for i, _fli in enumerate(fl):
            
            # get emissions data
            fl[i] = emi_model.eval(fl[i])

            # Iterate over the columns in the DataFrame
            for column in fl[i].dataframe.columns:
                # Replace NaN values in the column with the value from the previous row
                fl[i].dataframe[column] = fl[i].dataframe[column].fillna(method="ffill")

            # emission indices
            eis = {
                # primary combustion products
                "CO2": 3.16,
                "H2O": 1.23,
                "SO2": 0.00084,
                # secondary combustion products
                "nvPM": fl[i]["nvpm_ei_m"],
                "NO": 0.95 * fl[i]["nox_ei"],
                "NO2": 0.05 * fl[i]["nox_ei"],
                "CO": fl[i]["co_ei"],
                # hydrocarbon speciation
                "HCHO": 0.12 * fl[i]["hc_ei"],  # formaldehyde
                "CH3CHO": 0.04 * fl[i]["hc_ei"],  # acetaldehyde
                "C2H4": 0.15 * fl[i]["hc_ei"],  # ethylene
                "C3H6": 0.04 * fl[i]["hc_ei"],  # propene
                "C2H2": 0.04 * fl[i]["hc_ei"],  # acetylene
                "BENZENE": 0.02 * fl[i]["hc_ei"],  # benzene
            }

            # calculate emission mass total at each waypoint
            for species, ei in eis.items():
                fl[i][species] = ei * fl[i]["fuel_burn"] # [kg]

        return fl

    def sim_plumes(self) -> list[pd.DataFrame]:
        """Simulate plume dispersion/advection using Pycontrails Dry Advection Model."""
        plume_params = self.plume_params
        sim_params = self.sim_params
        met = self.met
        fl = self.fl

        dry_adv = DryAdvection(
            met,
            max_age=sim_params.t_pl[2],
            dt_integration=sim_params.t_pl[1],
            shear=plume_params.shear,
        )

        pl = []

        for i, fli in enumerate(fl):
            pli = dry_adv.eval(fli)
            pl.append(pli)

            # convert both flights and plumes to dataframes
            fl[i] = fl[i].dataframe
            for column in fl[i].columns:
                # Replace NaN values in the column with the value from the previous row
                fl[i][column] = fl[i][column].fillna(method="ffill")

            pl[i] = pl[i].dataframe

            # calc plume heading
            pl[i] = calc_heading(pl[i])
            pl[i]["flight_id"] = fl[i]["flight_id"][0]
            fl[i]["waypoint"] = fl[i].index

        # concatenate all flights and plumes into single dfs
        fl_df = pd.concat(fl)
        pl_df = pd.concat(pl)

        # merge the two dataframes
        fl = fl_df
        pl = pd.merge(
            fl_df[
                [
                    "flight_id",
                    "waypoint",
                    "fuel_flow",
                    "fuel_burn",
                    "true_airspeed",
                    "CO2",
                    "H2O",
                    "SO2",
                    "NO",
                    "NO2",
                    "CO",
                    "HCHO",
                    "CH3CHO",
                    "C2H4",
                    "C3H6",
                    "C2H2",
                    "BENZENE",
                    "nvPM",
                ]
            ],
            pl_df[
                [
                    "flight_id",
                    "waypoint",
                    "time",
                    "age",
                    "longitude",
                    "latitude",
                    "level",
                    "width",
                    "depth",
                    "heading",
                    "sigma_yy",
                    "sigma_yz",
                    "sigma_zz",
                ]
            ],
            on=["flight_id", "waypoint"],
        ).sort_values(by=["time", "flight_id", "waypoint"])

        pl["sin_a"] = np.sin(np.radians(pl["heading"]))
        pl["cos_a"] = np.cos(np.radians(pl["heading"]))
        pl["altitude"] = units.pl_to_m(pl["level"])
        pl["time"] = pl["time"] - sim_params.t_pl[1]

        return fl, pl

    def run_boxm(self) -> xr.Dataset:
        """Run BOXM."""
        # Initialize the box model dataset
        self.gen_inputs()
        self.gen_outputs()
        # Run the box model
        # self.run_boxm_f90()
        # # Unstack the box model dataset
        # self.unstack()

        # return self.boxm_ds_unstacked

    def gen_inputs(self):
        """Generate BOXM inputs."""
        # Initialize the box model dataset
        self.init_boxm_nc()
        # Initialise flight dataset
        self.init_fl_nc()
        #Initialise plume dataset
        self.init_pl_nc()
        

    def gen_outputs(self):
        """Generate BOXM output templates."""
        # Initialise boxm coarse output dataset
        self.init_boxm_out_nc()
        # Initialise patch table output dataset
        self.init_patch_table_nc()
        # Initialise plume output dataset
        self.init_pl_out_nc()  

    # Methods for running the box model

    def init_boxm_nc(self):
        """Initialize the box model dataset (met + background chem on coarse grid)."""

        # --- Merge meteorology and background chemistry fields ---
        self.boxm_ds = xr.merge([self.met.data, self.bg_chem])

        # --- Drop unneeded diagnostic fields ---
        drop_vars = [
            "specific_humidity",
            "relative_humidity",
            "eastward_wind",
            "northward_wind",
            "lagrangian_tendency_of_air_pressure",
            "month",
        ]
        self.boxm_ds = self.boxm_ds.drop_vars([v for v in drop_vars if v in self.boxm_ds.variables])

        # --- Assign useful metadata ---
        self.boxm_ds = self.boxm_ds.assign_attrs(
            ts_fl=self.sim_params.t_fl[1].total_seconds(),
            ts_pl=self.sim_params.t_pl[1].total_seconds(),
            ts_sim=self.sim_params.t_sim[1].total_seconds(),
            hres_sim_c=self.sim_params.hres_sim_c,
            vres_sim_c=self.sim_params.vres_sim_c,
            hres_sim_f=self.sim_params.hres_sim_f,
            vres_sim_f=self.sim_params.vres_sim_f,
            species_in=self.chem_params.species_in,
            species_out=self.chem_params.species_out,
            species_out_num=self.chem_params.species_out_num,

            description="BOXM coarse-grid meteorology and background chemistry fields",
            note="Emissions and plume segments handled separately via PL.NC and FL.NC",
        )

        # Flatten spatial dimensions for easy Fortran indexing
        # (Fortran expects a 1D cell index)
        self.boxm_ds_stacked = self.boxm_ds.stack(
            {"cell": ["level", "longitude", "latitude"]}
        ).reset_index("cell")

        # Delete any existing NetCDF
        nc_path = pathlib.Path(f"{self.inputs_job}/boxm.nc")
        if nc_path.exists():
            print("Deleting existing boxm.nc")
            nc_path.unlink()

        # Save to NetCDF file
        self.boxm_ds_stacked.to_netcdf(nc_path, mode="w")
        print(f"Saved {nc_path}")

    def init_fl_nc(self):
        """Initialize the flight dataset for BOXM."""
        df = self.fl.copy()

        # Set multi-index (flight_id, waypoint) and convert to xarray
        fl_ds = df.set_index(["flight_id", "waypoint"]).to_xarray()

        # Drop individual species variables
        fl_ds = fl_ds.drop_vars(["air_temperature", "specific_humidity", 
                                 "nox_ei", "co_ei", "hc_ei", "nvpm_ei_m",
                                 "nvpm_ei_n", "co2", "h2o", "so2", 
                                 "sulphates", "oc", "nox", "co", "hc", 
                                 "nvpm_mass", "nvpm_number", "nvPM"])

        # Save to NetCDF
        nc_path = pathlib.Path(f"{self.inputs_job}/fl.nc")
        if nc_path.exists():
            nc_path.unlink()
        fl_ds.to_netcdf(nc_path, mode="w")
        print(f"Saved {nc_path}")
    
    def init_pl_nc(self):
        """Initialize the plume dataset for the box model (PL.NC)."""
        df = self.pl.copy()

        species_cols = self.chem_params.species_in
        all_species_cols = ['CO2', 'H2O', 'SO2', 'NO', 
                            'NO2', 'CO', 'HCHO', 'CH3CHO', 
                            'C2H4', 'C3H6', 'C2H2', 'BENZENE', 'nvPM']
        # Set multi-index (flight_id, waypoint, time)
        pl_ds = df.set_index(["flight_id", "waypoint", "time"]).to_xarray()

        # Stack species into a single 4D DataArray (flight_id, waypoint, time, species)
        species_data = []
        for col in species_cols:
            species_data.append(pl_ds[col])
        
        pl_ds["emi_species_mass"] = xr.concat(
            species_data,
            dim=pd.Index(species_cols, name="species")
        ).transpose("flight_id", "waypoint", "time", "species")  # <-- ensure correct order
        
        # Drop ALL emission species variables (both stacked and unstacked)
        vars_to_drop = [col for col in all_species_cols if col in pl_ds.data_vars]
        if vars_to_drop:
            pl_ds = pl_ds.drop_vars(vars_to_drop)

        # drop remaining unneeded variables
        # pl_ds = pl_ds.drop_vars(["sin_a", "cos_a", "

        # Save to NetCDF
        nc_path = pathlib.Path(f"{self.inputs_job}/pl.nc")
        if nc_path.exists():
            nc_path.unlink()
        pl_ds.to_netcdf(nc_path, mode="w")
        print(f"Saved {nc_path}")


    def init_boxm_out_nc(self):
        """Initialize the box model coarse output dataset (BOXM_C_OUT.NC)."""
        species_out = np.array(self.chem_params.species_out, dtype="U10")

        self.boxm_out = xr.Dataset(
            data_vars={
                "Y_bg_c": (
                    ("time", "level", "longitude", "latitude", "species_out"),
                    da.zeros(
                        (
                            len(self.times_sim), 
                            len(self.levels), 
                            len(self.lons), 
                            len(self.lats),
                            len(species_out)
                        ),
                        dtype=float
                    ),
                    {"units": "mol_cm3"},
                ),
                "Y_del_c": (
                    ("time", "level", "longitude", "latitude", "species_out"),
                    da.zeros(
                        (
                            len(self.times_sim), 
                            len(self.levels), 
                            len(self.lons), 
                            len(self.lats),
                            len(species_out)
                        ),
                        dtype=float
                    ),
                    {"units": "mol_cm3"},
                ),
                "active_flag": (
                    ("time", "level", "longitude", "latitude"),
                    da.zeros(
                        (
                            len(self.times_sim), 
                            len(self.levels), 
                            len(self.lons), 
                            len(self.lats)
                        ),
                        dtype=bool
                    )
                ),
            },
            coords={
                "time": self.times_sim,
                "level": self.levels,
                "longitude": self.lons,
                "latitude": self.lats,
                "species_out": species_out,
            }
        )

        # Flatten spatial dimensions for easy Fortran indexing
        # (Fortran expects a 1D cell index)
        self.boxm_out_stacked = self.boxm_out.stack(
            {"cell": ["level", "longitude", "latitude"]}
        ).reset_index("cell")

        # Delete any existing NetCDF
        nc_path = pathlib.Path(f"{self.outputs_job}/boxm_out.nc")
        if nc_path.exists():
            print("Deleting existing boxm_out.nc")
            nc_path.unlink()

        # Save to NetCDF file
        self.boxm_out_stacked.to_netcdf(nc_path, mode="w")
        print(f"Saved {nc_path}")
        
    def init_patch_table_nc(self):
        """Initialize the patch output dataset (PATCH_OUT.NC)."""
        species_out = np.array(self.chem_params.species_out, dtype="U10")

        self.patch_table = xr.Dataset(
            data_vars={
                "Y_del_f": (
                    ("row", "species_out"), 
                    da.zeros(
                        (
                            0, 
                            len(species_out)
                        ),
                        dtype=float,
                    ),  
                    {"units": "mol_cm3"},
                ),
            },
            coords={
            "row": ("row", np.array([], dtype=int)),   # 0-length row dimension
            "time": ("row", np.array([], dtype="datetime64[ns]")),
            "patch_id": ("row", np.array([], dtype=int)),

            "latitude": ("row", np.array([], dtype=float)),
            "longitude": ("row", np.array([], dtype=float)),
            "level": ("row", np.array([], dtype=int)),

            "latitude_f": ("row", np.array([], dtype=float)),
            "longitude_f": ("row", np.array([], dtype=float)),
            "level_f": ("row", np.array([], dtype=int)),

            "species_out": ("species_out", species_out),
            }
        )

        # Delete any existing NetCDF
        nc_path = pathlib.Path(f"{self.outputs_job}/patch_table.nc")
        if nc_path.exists():
            print("Deleting existing patch_table.nc")
            nc_path.unlink()

        # Save to NetCDF file
        self.patch_table.to_netcdf(nc_path, mode="w")
        print(f"Saved {nc_path}")
        
    def init_pl_out_nc(self):
        """Initialize the plume output dataset (PL_OUT.NC)."""
        # Load the input plume dataset
        pl_ds = xr.open_dataset(f"{self.inputs_job}/pl.nc")
        
        # Create output dataset with same structure
        self.pl_out = pl_ds.copy(deep=True)
        
        # Rename and repurpose the species mass variable
        # Store DELTA mass (change from initial emissions due to chemistry)
        self.pl_out = self.pl_out.rename({"emi_species_mass": "delta_species_mass"})
        
        # Zero out delta mass (initially, no chemistry has occurred)
        self.pl_out["delta_species_mass"].values[:] = 0.0
        
        # If output species differ from input species, create new variable
        if set(self.chem_params.species_in) != set(self.chem_params.species_out):
            self.pl_out = self.pl_out.drop_vars("species")
            
            species_out = np.array(self.chem_params.species_out, dtype="U10")
            
            # Create delta_species_mass for output species
            self.pl_out["delta_species_mass"] = (
                ("flight_id", "waypoint", "time", "species_out"),
                np.zeros((
                    len(self.pl_out.flight_id),
                    len(self.pl_out.waypoint),
                    len(self.pl_out.time),
                    len(species_out)
                ), dtype=float),
                {
                    "units": "kg",
                    "long_name": "Change in species mass due to chemistry",
                    "note": "Add to initial emission mass (from fl.nc) to get total mass"
                }
            )
            self.pl_out["species_out"] = species_out
        else:
            # Just update metadata
            self.pl_out["delta_species_mass"].attrs = {
                "units": "kg",
                "long_name": "Change in species mass due to chemistry",
                "note": "Add to initial emission mass (from fl.nc) to get total mass"
            }

        # Add aspect ratio column
        self.pl_out["aspect_ratio"] = (
            self.pl_out["width"] / self.pl_out["depth"]
        )
        
        # Add metadata
        self.pl_out.attrs.update({
            "description": "BOXM plume chemistry output (mass deltas)",
            "created": pd.Timestamp.now().isoformat(),
            "note": "Geometry copied from pl.nc. Total mass = emi_species_mass (fl.nc) + delta_species_mass (pl_out.nc)"
        })
        
        # Save
        nc_path = pathlib.Path(f"{self.outputs_job}/pl_out.nc")
        if nc_path.exists():
            nc_path.unlink()
        self.pl_out.to_netcdf(nc_path, mode="w")
        print(f"Saved {nc_path}")

    def run_boxm_f90(self):
        """Run the box model in fortran using subprocess."""
        # Run the box model
        subprocess.call(
            [self.run_path + "boxm", self.job_id],
        )

        # open nc file
        # self.boxm_ds = xr.open_dataset(f"{self.inputs_job}/boxm.nc")

    def unstack(self):
        """Unstack the box model dataset."""
        print("Chunking the dataset")
        # Convert the dataset to a Dask dataset
        self.boxm_ds = self.boxm_ds.chunk({"cell": 100})  # Adjust chunk size based on your mem

        print("Set coords")
        # Convert 'level', 'lat', and 'lon' to coordinates
        self.boxm_ds_unstacked = self.boxm_ds.set_coords(["level", "longitude", "latitude"])

        print("Set index")
        # Create a multi-index for the 'cell' dimension
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.set_index(
            cell=["level", "longitude", "latitude"]
        )

        print("Unstack the dataset")
        # Unstack the dataset
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.unstack("cell")

        print("Compute the result")
        # # Compute the result to trigger the lazy evaluation
        # self.boxm_ds_unstacked = self.boxm_ds_unstacked.compute()


# Functions used in GPAT Model
def grab_species_num(run_path, species_out: np.array) -> np.array:
    """Grab the species numbers for the species of interest in output."""
    # Read species names from the file into a list
    with open(f"{run_path}species_num.txt") as file:
        species_list = [line.strip() for line in file]

    # Create a dictionary mapping species names to their line numbers
    species_dict = {species: index for index, species in enumerate(species_list)}

    return np.array([species_dict[species] + 1 for species in species_out])


def calc_heading(pl_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate heading for each plume.

    Parameters
    ----------
    pl_df : pd.DataFrame
        DataFrame containing plume data.

    Returns
    -------
    pd.DataFrame
        DataFrame containing plume data with heading.
    """
    # Sort the dataframe by time and waypoint
    pl_df = pl_df.sort_values(by=["time", "waypoint"])

    # Group the dataframe by the timestep and apply the function
    pl_df["heading"] = pl_df.groupby("time").apply(calculate_heading_g).reset_index(drop=True)

    return pl_df


def calculate_heading_g(group):
    """Calculate heading for each timestep.

    Parameters
    ----------
    group : pd.DataFrame
        DataFrame containing plume data for a single timestep.

    Returns
    -------
    pd.Series
        Series containing heading for each plume in the timestep.
    """

    g = Geod(ellps="WGS84")

    startlat = group["latitude"].values[:-1]
    startlon = group["longitude"].values[:-1]
    endlat = group["latitude"].values[1:]
    endlon = group["longitude"].values[1:]
    az12, az21, dist = g.inv(startlon, startlat, endlon, endlat)

    heading = (90 - az12) % 360

    return pd.Series(
        np.concatenate([[heading[0]], heading]) if len(heading) > 0 else [np.nan], index=group.index
    )


def calc_continuous(plume: GeoVectorDataset):
    """Calculate the continuous segments of this timestep.

    Mutates parameter ``contrail`` in place by setting or updating the
    "continuous" variable.

    Parameters
    ----------
    contrail : GeoVectorDataset
        GeoVectorDataset instance onto which "continuous" is set.

    Raises
    ------
    ValueError
        If ``contrail`` is empty.
    """

    if not plume:
        raise ValueError("Cannot calculate continuous on an empty contrail")
    same_flight = plume["flight_id"][:-1] == plume["flight_id"][1:]
    consecutive_waypoint = np.diff(plume["waypoint"]) == 1
    continuous = np.empty(plume.size, dtype=bool)
    continuous[:-1] = same_flight & consecutive_waypoint
    continuous[-1] = False  # This fails if contrail is empty
    plume.update(continuous=continuous)  # overwrite continuous


def calc_sza(latitudes, longitudes, timesteps):
    """Calculate szas for each cell at all timesteps.

    Parameters
    ----------
    latitudes : np.array
        Array of latitudes.
    longitudes : np.array
        Array of longitudes.
    timesteps : np.array
        Array of timesteps.

    Returns
    -------
    np.array
        Array of szas for each cell at all timesteps.
    """

    sza = np.zeros((len(latitudes), len(longitudes), len(timesteps)))

    for lon, lonval in enumerate(longitudes):
        for lat, latval in enumerate(latitudes):
            theta_rad = geo.orbital_position(timesteps)

            sza[lat, lon, :] = np.arccos(
                geo.cosine_solar_zenith_angle(lonval, latval, timesteps, theta_rad)
            )
    return sza


# Validation
def mc_test(params, fl_df, pl_df, chem_ds):
    """Check if mass is conserved in the box model.

    Parameters
    ----------
    params : pd.Series
        Series containing simulation parameters.
    fl_df : pd.DataFrame
        DataFrame containing flight data.
    pl_df : pd.DataFrame
        DataFrame containing plume data.
    chem_ds : xr.Dataset
        Dataset containing chemical data.

    Returns
    -------
    pd.DataFrame
        DataFrame containing mass conservation data.
    """

    # Initialize the dictionary
    vecmass = {emi_species: [] for emi_species in chem_ds["emi_species"].values.tolist()}
    gridmass = {emi_species: [] for emi_species in chem_ds["emi_species"].values.tolist()}
    mc = {emi_species: [] for emi_species in chem_ds["emi_species"].values.tolist()}

    # Constants
    mm = [30.01, 46.01, 28.01, 30.03, 44.05, 28.05, 42.08, 26.04, 78.11]  # g/mol
    NA = 6.022e23  # Avogadro's number

    for s, emi_species in enumerate(chem_ds["emi_species"].values):
        max_fl_time = fl_df["time"].max()

        for ts, t in enumerate(pl_df["time"].unique()[:-1]):
            if ts == 0:
                total_vector_mass = 0
                total_grid_mass = 0
                percent_mass_conserved = 0
                vecmass[emi_species].append(total_vector_mass)
                gridmass[emi_species].append(total_grid_mass)
                mc[emi_species].append(percent_mass_conserved)
                continue

            previous_time = pl_df["time"].unique()[ts - 1]
            fl_snapshot = fl_df[fl_df["time"] == previous_time]

            if t <= max_fl_time:
                # Accumulate vector mass for all flights
                vector_mass = fl_snapshot[emi_species]

                total_vector_mass += vector_mass.sum()

            # Grab plume mass from grid data
            grid_concs = chem_ds["emi"].sel(emi_species=emi_species, time=t)

            if (grid_concs == 0).all():
                pass
            else:
                # Compute the boolean indexer first
                grid_concs_over_zero = grid_concs > 0
                grid_concs_over_zero = grid_concs.where(grid_concs_over_zero, drop=True)

                grid_mass = (
                    grid_concs_over_zero
                    * chem_ds["M"].sel(time=t)
                    * 1e-9
                    * (mm[s] / NA)
                    * params.loc["vres_sim"]
                    * units.latitude_distance_to_m(params.loc["hres_sim"])
                    * units.longitude_distance_to_m(
                        params.loc["hres_sim"],
                        (params.loc["lat_bounds"][0] + params.loc["lat_bounds"][1]) / 2,
                    )
                    * 1e03
                )  # convert to kg/m^3

                total_grid_mass = grid_mass.sum().item()

                percent_mass_conserved = total_grid_mass / total_vector_mass * 100

            # Append the percentage to the list in the dictionary
            vecmass[emi_species].append(total_vector_mass)
            gridmass[emi_species].append(total_grid_mass)
            mc[emi_species].append(percent_mass_conserved)

    # convert the dictionary to a DataFrame
    vecmass = pd.DataFrame(
        vecmass, index=pl_df["time"].unique()[:-1], columns=chem_ds["emi_species"].values.tolist()
    )
    gridmass = pd.DataFrame(
        gridmass, index=pl_df["time"].unique()[:-1], columns=chem_ds["emi_species"].values.tolist()
    )
    mc = pd.DataFrame(
        mc, index=pl_df["time"].unique()[:-1], columns=chem_ds["emi_species"].values.tolist()
    )

    return vecmass, gridmass, mc


def boxm_test(run_path, data_path, job_id, chem_ds):
    """Run the box model for selected cells and job_id.

    Parameters
    ----------
    run_path : str
        Path to the box model executable.
    data_path : str
        Path to the data directory.
    job_id : str
        Job ID.
    chem_ds : xr.Dataset
        Dataset containing chemical data.

    Returns
    -------
    xr.Dataset
        Dataset containing chemical data for the selected cells.
    """

    # create input file for original boxm
    gen_boxm_orig_input(data_path, chem_ds, job_id)

    gen_zen_file(data_path, chem_ds, job_id)

    gen_emi_file(data_path, chem_ds, job_id)

    # # calls fortran with input file and generates .OUT files
    subprocess.call(
        [run_path + "boxm_orig", data_path, job_id],
    )

    return update_chem_ds(data_path, chem_ds, job_id)


def gen_boxm_orig_input(data_path, cell_chem_ds, job_id):
    """Generate the input file for the original box model.

    Parameters
    ----------
    data_path : str
        Path to the data directory.
    cell_chem_ds : xr.Dataset
        Dataset containing chemical data.
    job_id : str
        Job ID.
    """

    # delete any existing input files
    if pathlib.Path(f"inputs/{job_id}/boxm_input.txt").exists():
        pathlib.Path(f"inputs/{job_id}/boxm_input.txt").unlink()

    # open file using a context manager
    with open(f"inputs/{job_id}/boxm_input.txt", "w") as boxm_input:
        start_time = pd.to_datetime(cell_chem_ds["time"].values[0])
        end_time = pd.to_datetime(cell_chem_ds["time"].values[-1])
        runtime = int((end_time - start_time) / np.timedelta64(1, "D")) % 365
        day = start_time.day
        month = start_time.month
        year = start_time.year
        altitude = cell_chem_ds["altitude"].item()
        plevel = cell_chem_ds["level"].item()
        level = get_pressure_level(altitude)
        longitude = cell_chem_ds["longitude"].item()
        longbox = longitude_to_longbox(longitude)
        latitude = cell_chem_ds["latitude"].item()
        latbox = latitude_to_latbox(latitude)
        M = cell_chem_ds["M"].values[0]
        # P = cell_chem_ds["air_pressure"].item()
        H2O = cell_chem_ds["H2O"].values[0]
        temp = cell_chem_ds["air_temperature"].values[0]

        boxm_input.write(
            f"{day}\n{month}\n{year}\n{level}\n{longbox}\n{latbox}\n{runtime}\n{M}\n{plevel}"
            f"\n{H2O}\n{temp}\n"
        )
    for s in [
        "NO2",
        "NO",
        "O3",
        "CO",
        "CH4",
        "HCHO",
        "CH3CHO",
        "CH3COCH3",
        "C2H6",
        "C2H4",
        "C3H8",
        "C3H6",
        "C2H2",
        "NC4H10",
        "TBUT2ENE",
        "BENZENE",
        "TOLUENE",
        "OXYL",
        "C5H8",
        "H2O2",
        "HNO3",
        "C2H5CHO",
        "CH3OH",
        "MEK",
        "CH3OOH",
        "PAN",
        "MPAN",
    ]:
        boxm_input.write(f"{cell_chem_ds['bg_chem'].sel(species=s).item()}\n")

    boxm_input.close()


def gen_zen_file(data_path, cell_chem_ds, job_id):
    """Generate the ZEN file for the original box model.

    Parameters
    ----------
    data_path : str
        Path to the data directory.
    cell_chem_ds : xr.Dataset
        Dataset containing chemical data.
    job_id : str
        Job ID.
    """

    # delete any existing input files
    zen_file_path = pathlib.Path(f"{data_path}outputs/{job_id}/zen.csv")
    if zen_file_path.exists():
        zen_file_path.unlink()

    # Extract the sza data and convert it to a DataFrame
    sza_data = cell_chem_ds["sza"].values
    sza_df = pd.DataFrame(sza_data, columns=["sza"])

    # Write the DataFrame to a CSV file
    sza_df.to_csv(zen_file_path, index=False, header=False)


def gen_emi_file(data_path, cell_chem_ds, job_id):
    """Generate the EMI file for the original box model.

    Parameters
    ----------
    data_path : str
        Path to the data directory.
    cell_chem_ds : xr.Dataset
        Dataset containing chemical data.
    job_id : str
        Job ID.
    """

    # delete any existing input files
    emi_file_path = pathlib.Path(f"{data_path}outputs/{job_id}/emi.csv")
    if emi_file_path.exists():
        emi_file_path.unlink()

    # Extract the emi data and convert it to a DataFrame
    emi_data = cell_chem_ds["emi"].values
    emi_df = pd.DataFrame(emi_data, columns=cell_chem_ds["emi_species"].values)

    # Write the DataFrame to a CSV file
    emi_df.to_csv(emi_file_path, index=False, header=False)


def latitude_to_latbox(latitude):
    """Convert latitude to latbox.

    Parameters
    ----------
    latitude : float
        Latitude.

    Returns
    -------
    int
        Latbox.
    """
    # Map the latitude to the range 0-1
    normalized_latitude = (latitude + 87.5) / 180

    # Map the normalized latitude to the range 1-72
    latbox = normalized_latitude * 36 + 1

    # Round to the nearest integer and return
    return round(latbox)


def longitude_to_longbox(longitude):
    """Convert longitude to longbox.

    Parameters
    ----------
    longitude : float
        Longitude.

    Returns
    -------
    int
        Longbox.
    """
    # Map the longitude to the range 0-1
    normalized_longitude = (longitude + 177.5) / 360

    # Map the normalized longitude to the range 1-144
    longbox = normalized_longitude * 72 + 1

    # Round to the nearest integer and return
    return round(longbox)


def get_pressure_level(alt):
    """Get the pressure level for a given altitude.

    Parameters
    ----------
    alt : float
        Altitude.

    Returns
    -------
    int
        Pressure level.
    """

    # Convert alt to pressure level (hPa)``
    chem_pressure_levels = np.array([962, 861, 759, 658, 556, 454, 353, 251, 150.5])

    # Convert altitude to pressure using a standard atmosphere model
    pressure = units.m_to_pl(alt)

    # Find the index of the closest value in the array
    return (np.abs(chem_pressure_levels - pressure)).argmin()


def update_chem_ds(data_path, cell_chem_ds, job_id):
    """Update the chemical dataset with the box model output.

    Parameters
    ----------
    data_path : str
        Path to the data directory.
    cell_chem_ds : xr.Dataset
        Dataset containing chemical data.
    job_id : str
        Job ID.

    Returns
    -------
    xr.Dataset
        Dataset containing updated chemical data.
    """

    sza_df = pd.read_csv(
        f"{data_path}outputs/{job_id}/ZEN.OUT", header=0, names=["TIME", "ZEN"], dtype=np.float64
    )

    J_df = pd.read_csv(
        f"{data_path}outputs/{job_id}/J.OUT",
        header=0,
        names=[
            "TIME",
            "J1",
            "J2",
            "J3",
            "J4",
            "J5",
            "J6",
            "J7",
            "J8",
            "J9",
            "J10",
            "J11",
            "J12",
            "J13",
            "J14",
            "J15",
            "J16",
            "J17",
            "J18",
            "J19",
            "J20",
            "J21",
            "J22",
            "J23",
            "J24",
            "J25",
            "J26",
            "J27",
            "J28",
            "J29",
            "J30",
            "J31",
            "J32",
            "J33",
            "J34",
            "J35",
            "J36",
            "J37",
            "J38",
            "J39",
            "J40",
            "J41",
            "J42",
            "J43",
            "J44",
            "J45",
            "J46",
            "J47",
            "J48",
            "J49",
            "J50",
        ],
        dtype=np.float64,
    )

    DJ_df = pd.read_csv(
        f"{data_path}outputs/{job_id}/DJ.OUT",
        header=0,
        names=[
            "TIME",
            "DJ1",
            "DJ2",
            "DJ3",
            "DJ4",
            "DJ5",
            "DJ6",
            "DJ7",
            "DJ8",
            "DJ9",
            "DJ10",
            "DJ11",
            "DJ12",
            "DJ13",
            "DJ14",
            "DJ15",
            "DJ16",
            "DJ17",
            "DJ18",
            "DJ19",
            "DJ20",
            "DJ21",
            "DJ22",
            "DJ23",
            "DJ24",
            "DJ25",
            "DJ26",
            "DJ27",
            "DJ28",
            "DJ29",
            "DJ30",
            "DJ31",
            "DJ32",
            "DJ33",
            "DJ34",
            "DJ35",
            "DJ36",
            "DJ37",
            "DJ38",
            "DJ39",
            "DJ40",
            "DJ41",
            "DJ42",
            "DJ43",
            "DJ44",
            "DJ45",
            "DJ46",
            "DJ47",
            "DJ48",
            "DJ49",
            "DJ50",
        ],
        dtype=np.float64,
    )

    RC_df = pd.read_csv(
        f"{data_path}outputs/{job_id}/RC.OUT",
        header=0,
        names=[
            "TIME",
            "RC1",
            "RC2",
            "RC3",
            "RC4",
            "RC5",
            "RC6",
            "RC7",
            "RC8",
            "RC9",
            "RC10",
            "RC11",
            "RC12",
            "RC13",
            "RC14",
            "RC15",
            "RC16",
            "RC17",
            "RC18",
            "RC19",
            "RC20",
            "RC21",
            "RC22",
            "RC23",
            "RC24",
            "RC25",
            "RC26",
            "RC27",
            "RC28",
            "RC29",
            "RC30",
            "RC31",
            "RC32",
            "RC33",
            "RC34",
            "RC35",
            "RC36",
            "RC37",
            "RC38",
            "RC39",
            "RC40",
            "RC41",
            "RC42",
            "RC43",
            "RC44",
            "RC45",
            "RC46",
            "RC47",
            "RC48",
            "RC49",
            "RC50",
        ],
        dtype=np.float64,
    )

    # get species names
    header_names = ["TIME", *list(cell_chem_ds["species"].values)]

    Y_df = pd.read_csv(
        f"{data_path}/outputs/{job_id}/Y.OUT", header=0, names=header_names, dtype=np.float64
    )

    # # Update the chem_ds_stacked with the new data
    # Update zen data
    cell_chem_ds["sza_orig"] = (["time"], da.zeros(cell_chem_ds.sizes["time"]))
    cell_chem_ds["sza_orig"].loc[:] = sza_df["ZEN"].values * np.pi / 180

    cell_chem_ds["J_orig"] = (["time", "photol_params"], da.zeros((cell_chem_ds.sizes["time"], 5)))
    for pp, photol_params in enumerate(J_df.columns[1:6]):
        cell_chem_ds["J_orig"].loc[:, pp] = J_df[photol_params].values

    cell_chem_ds["DJ_orig"] = (["time", "photol_coeffs"], da.zeros((cell_chem_ds.sizes["time"], 5)))
    for pc, photol_coeffs in enumerate(DJ_df.columns[1:6]):
        cell_chem_ds["DJ_orig"].loc[:, pc] = DJ_df[photol_coeffs].values

    cell_chem_ds["RC_orig"] = (["time", "therm_coeffs"], da.zeros((cell_chem_ds.sizes["time"], 5)))
    for tc, therm_coeffs in enumerate(RC_df.columns[1:6]):
        cell_chem_ds["RC_orig"].loc[:, tc] = RC_df[therm_coeffs].values

    cell_chem_ds["Y_orig"] = (
        ["time", "species_out"],
        da.zeros((cell_chem_ds.sizes["time"], cell_chem_ds.sizes["species_out"])),
    )
    for _s, species_out in enumerate(cell_chem_ds["species_out"].values):
        cell_chem_ds["Y_orig"].loc[:, species_out] = Y_df[species_out].values

    return cell_chem_ds


def parse_args():
    """Parse command line arguments.

    Returns
    -------
    argparse.Namespace
        Namespace containing command line arguments.
    """

    parser = argparse.ArgumentParser(description="Overwrite parameters from command line")

    # FlParams arguments
    parser.add_argument("--t0_fl", type=str, help="Flight start time")
    parser.add_argument("--rt_fl", type=int, help="Flight run time in minutes")
    parser.add_argument("--ts_fl", type=int, help="Flight time step in minutes")
    parser.add_argument("--ac_type", type=str, help="Aircraft type")
    parser.add_argument("--fl0_speed", type=float, help="Flight speed in m/s")
    parser.add_argument("--fl0_heading", type=float, help="Flight heading in degrees")
    parser.add_argument("--fl0_coords0", type=str, help="Flight coordinates (lat, lon, alt)")
    parser.add_argument("--sep_dist", type=str, help="Separation distance (dx, dy, dz)")
    parser.add_argument("--n_ac", type=int, help="Number of aircraft")

    # PlumeParams arguments
    parser.add_argument("--dt_integration", type=int, help="Integration time step in minutes")
    parser.add_argument("--max_age", type=str, help="Maximum age of the plume in hours")
    parser.add_argument("--depth", type=float, help="Initial plume depth in meters")
    parser.add_argument("--width", type=float, help="Initial plume width in meters")
    parser.add_argument(
        "--hres_pl", type=float, help="Horizontal resolution of the plume in degrees"
    )
    parser.add_argument("--vres_pl", type=float, help="Vertical resolution of the plume in meters")
    parser.add_argument("--n_slices", type=int, help="Number of slices")

    # SimParams arguments
    parser.add_argument("--t0_sim", type=str, help="Simulation start time")
    parser.add_argument("--rt_sim", type=int, help="Simulation runtime in hours")
    parser.add_argument("--ts_sim", type=int, help="Simulation time step in seconds")
    parser.add_argument("--lat_bounds", type=str, help="Latitude bounds (min, max)")
    parser.add_argument("--lon_bounds", type=str, help="Longitude bounds (min, max)")
    parser.add_argument("--alt_bounds", type=str, help="Altitude bounds (min, max)")
    parser.add_argument("--hres_sim", type=float, help="Horizontal resolution in degrees")
    parser.add_argument("--vres_sim", type=float, help="Vertical resolution in meters")
    parser.add_argument("--eastward_wind", type=float, help="Eastward wind in m/s")
    parser.add_argument("--northward_wind", type=float, help="Northward wind in m/s")
    parser.add_argument(
        "--lagrangian_tendency_of_air_pressure",
        type=float,
        help="Lagrangian tendency of air pressure in m/s",
    )
    parser.add_argument("--species_in", type=str, help="Input species (comma-separated)")
    parser.add_argument("--species_out", type=str, help="Output species (comma-separated)")
    parser.add_argument("--gpat_path", type=str, help="GPAT directory")
    parser.add_argument("--job_id", type=str, help="Job ID")
    parser.add_argument("--run_gpat", action="store_true", help="Run the GPAT model")

    return parser.parse_args()


def update_fl_params_from_args(params, args):
    """Update FlParams from command line arguments.

    Parameters
    ----------
    params : FlParams
        FlParams instance.
    args : argparse.Namespace
        Namespace containing command line arguments.
    """
    if args.t0_fl:
        params.t0_fl = pd.to_datetime(args.t0_fl)
    if args.rt_fl:
        params.rt_fl = pd.Timedelta(minutes=args.rt_fl)
    if args.ts_fl:
        params.ts_fl = pd.Timedelta(minutes=args.ts_fl)
    if args.ac_type:
        params.ac_type = args.ac_type
    if args.fl0_speed:
        params.fl0_speed = args.fl0_speed
    if args.fl0_heading:
        params.fl0_heading = args.fl0_heading
    if args.fl0_coords0:
        params.fl0_coords0 = tuple(map(float, args.fl0_coords0.split(",")))
    if args.sep_dist:
        params.sep_dist = tuple(map(float, args.sep_dist.split(",")))
    if args.n_ac:
        params.n_ac = args.n_ac


def update_plume_params_from_args(params, args):
    """Update PlParams from command line arguments.

    Parameters
    ----------
    params : PlParams
        PlParams instance.
    args : argparse.Namespace
        Namespace containing command line arguments.
    """
    if args.dt_integration:
        params.dt_integration = pd.Timedelta(minutes=args.dt_integration)
    if args.max_age == "ID":
        params.max_age = args.max_age
    elif args.max_age:
        params.max_age = pd.Timedelta(hours=int(args.max_age))
    if args.depth:
        params.depth = args.depth
    if args.width:
        params.width = args.width
    if args.hres_pl:
        params.hres_pl = args.hres_pl
    if args.vres_pl:
        params.vres_pl = args.vres_pl
    if args.n_slices:
        params.n_slices = args.n_slices


def update_sim_params_from_args(params, args):
    """Update SimParams from command line arguments.

    Parameters
    ----------
    params : SimParams
        SimParams instance.
    args : argparse.Namespace
        Namespace containing command line arguments.
    """

    if args.t0_sim:
        params.t0_sim = pd.to_datetime(args.t0_sim)
    if args.rt_sim:
        params.rt_sim = pd.Timedelta(hours=args.rt_sim)
    if args.ts_sim:
        params.ts_sim = pd.Timedelta(seconds=args.ts_sim)
    if args.lat_bounds:
        params.lat_bounds = tuple(map(float, args.lat_bounds.split(",")))
    if args.lon_bounds:
        params.lon_bounds = tuple(map(float, args.lon_bounds.split(",")))
    if args.alt_bounds:
        params.alt_bounds = tuple(map(float, args.alt_bounds.split(",")))
    if args.hres_sim:
        params.hres_sim = args.hres_sim
    if args.vres_sim:
        params.vres_sim = args.vres_sim
    if args.eastward_wind:
        params.eastward_wind = args.eastward_wind
    if args.northward_wind:
        params.northward_wind = args.northward_wind
    if args.lagrangian_tendency_of_air_pressure:
        params.lagrangian_tendency_of_air_pressure = args.lagrangian_tendency_of_air_pressure
    if args.species_in:
        params.species_in = tuple(args.species_in.split(","))
    if args.species_out:
        params.species_out = tuple(args.species_out.split(","))
    if args.gpat_path:
        params.gpat_path = args.gpat_path
    if args.job_id:
        params.job_id = args.job_id
    if args.run_gpat:
        params.run_gpat = args.run_gpat


# Function to convert dictionary to dataclass instance
def dict_to_dataclass(cls, dict_obj):
    """Convert a dictionary to a dataclass instance.

    Parameters
    ----------
    cls : dataclass
        Dataclass.
    dict_obj : dict
        Dictionary.

    Returns
    -------
    dataclass
        Dataclass instance.
    """

    return cls(**dict_obj)


# Function to filter out inherited params from ModelParams
def filter_inherited_params(instance, base_class):
    """Filter out inherited parameters from a dataclass instance.

    Parameters
    ----------
    instance : dataclass
        Dataclass instance.
    base_class : dataclass
        Base class.

    Returns
    -------
    dict
        Dictionary containing only the parameters unique to the instance.
    """

    base_fields = {f.name for f in fields(base_class)}
    instance_dict = asdict(instance)
    return {k: v for k, v in instance_dict.items() if k not in base_fields}
