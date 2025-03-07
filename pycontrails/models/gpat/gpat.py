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

import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Geod

from pycontrails.core import Flight, GeoVectorDataset, MetDataset, models
from pycontrails.core.models import Model
from pycontrails.models.dry_advection import DryAdvection
from pycontrails.models.emissions import Emissions
from pycontrails.models.gpat.plume_to_grid import plume_to_grid
from pycontrails.models.ps_model import PSFlight
from pycontrails.physics import constants, geo, thermo, units


### GPAT Model Parameters ###
@dataclass
class FlParams:
    """Default flight/fleet parameters."""

    t0_fl: pd.Timestamp = field(
        default_factory=lambda: pd.to_datetime("2022-01-20 13:00:00")
        )  # flight start time
    rt_fl: pd.Timedelta = field(
        default_factory=lambda: pd.Timedelta(minutes=60)
        )  # flight run time
    ts_fl: pd.Timedelta = field(
        default_factory=lambda: pd.Timedelta(minutes=2)
                                )  # flight time step
    ac_type: str = "A320"  # aircraft type
    fl0_speed: float = 100.0  # m/s
    fl0_heading: float = 0.0  # deg
    fl0_coords0: tuple[float, float, float] = (0.1, 0.125, 12500)  # lat, lon, alt [deg, deg, m]
    sep_dist: tuple[float, float, float] = (5000, 2000, 0)  # dx, dy, dz [m]
    n_ac: int = 1  # number of aircraft

@dataclass
class PlumeParams:
    """Default plume dispersion parameters."""

    dt_integration: pd.Timedelta = field(
        default_factory=lambda: pd.Timedelta(minutes=2)
        )  # integration time step
    max_age: str | pd.Timedelta = field(
        default_factory=lambda: pd.Timedelta(hours=2)
        )  # maximum age of the plume
    depth: float = 50.0  # initial plume depth, [m]
    width: float = 50.0  # initial plume width, [m]
    shear: float = 0.01  # wind shear [1/s]
    hres_pl: float = 0.01  # horizontal resolution of the plume, [deg]
    vres_pl: float = 500  # vertical resolution of the plume [m]
    n_slices: int = 10  # number of slices

@dataclass
class SimParams:
    """Default simulation parameters."""

    t0_sim: pd.Timestamp = field(
        default_factory=lambda: pd.to_datetime("2022-01-20 12:00:00")
        )  # simulation start time
    rt_sim: pd.Timedelta = field(
        default_factory=lambda: pd.Timedelta(hours=120)
        )  # simulation runtime
    ts_sim: pd.Timedelta = field(
        default_factory=lambda: pd.Timedelta(seconds=20)
        )  # simulation time step
    lat_bounds: tuple[float, float] = (0.0, 1.0)  # lat bounds [deg]
    lon_bounds: tuple[float, float] = (0.0, 1.0)  # lon bounds [deg]
    alt_bounds: tuple[float, float] = (12000, 13000)  # alt bounds [m]
    hres_sim: float = 0.01  # horizontal resolution [deg]
    vres_sim: float = 500  # vertical resolution [m]
    eastward_wind: float = 0.0  # m/s
    northward_wind: float = 0.0  # m/s
    lagrangian_tendency_of_air_pressure: float = 0.0  # m/s
    species_in: tuple = ("NO")
    species_out: tuple = ("O3", "NO2", "NO",
                          "NO3", "N2O5", "HNO3",
                          "HONO", "HO2", "OH",
                          "H2O2", "H2O", "CO",
                          "CH4", "CH3O2")
    run_path: str = None   
    data_path: str = None
    job_id: str = None
    run_gpat: bool = False
    date_created: pd.Timestamp = None
    species_out_num: tuple = None

class GPAT(Model):
    """Gridded Plume Analysis Tool (GPAT).

    Simulate aircraft trajectories, estimate aircraft performance, fuel burn and emissions. Then 
    aggregates emissions, bg chemistry and meteorology to an Eulerian grid for photochemical and 
    microphysical processing.

    Parameters
    ----------
    fl_params : FlParams
        Flight parameters.
    plume_params : PlumeParams
        Plume dispersion parameters.
    sim_params : SimParams
        Simulation parameters.
    """

    name = "GPAT"
    long_name = "Gridded Plume Analysis Tool"
    #default_params = (FlParams, PlumeParams, SimParams)

    def __init__(
            self,
            fl_params: FlParams,
            plume_params: PlumeParams,
            sim_params: SimParams
            ):
        super().__init__()

        # Generate the grid
        self.lats_pl = np.arange(
            sim_params.lat_bounds[0], 
            sim_params.lat_bounds[1] + plume_params.hres_pl / 2, 
            plume_params.hres_pl
        )
        self.lons_pl = np.arange(
            sim_params.lon_bounds[0], 
            sim_params.lon_bounds[1] + plume_params.hres_pl / 2, 
            plume_params.hres_pl
        )
        self.lats = np.arange(
            sim_params.lat_bounds[0], 
            sim_params.lat_bounds[1] + sim_params.hres_sim / 2, 
            sim_params.hres_sim
        )
        self.lons = np.arange(
            sim_params.lon_bounds[0], 
            sim_params.lon_bounds[1] + sim_params.hres_sim / 2, 
            sim_params.hres_sim
        )
        self.alts = np.arange(
            sim_params.alt_bounds[0], 
            sim_params.alt_bounds[1] + sim_params.vres_sim / 2, 
            sim_params.vres_sim
        )
        self.levels = units.m_to_pl(self.alts)

        self.times = pd.date_range(
            start=sim_params.t0_sim,
            end=sim_params.t0_sim + sim_params.rt_sim,
            freq=sim_params.ts_sim,
        )

        self.total_volume = (
            (self.lats_pl[-1] - self.lats_pl[0])
            * (self.lons_pl[-1] - self.lons_pl[0])
            * (self.alts[-1] - self.alts[0])
        )

        if plume_params.max_age == "ID":
                plume_params.max_age = plume_params.dt_integration

        if sim_params.run_path is None:
            self.run_path = os.environ['PYCONTRAILSDIR'] + "models/gpat/"

        else:
            self.run_path = sim_params.run_path

        if sim_params.data_path is None:
            self.data_path = "/projects/Impact_of_aviation_on_climate/Kieran2024/"

        else:
            self.data_path = sim_params.data_path

        if sim_params.job_id is None:
            try:
                self.job_id = os.environ['SLURM_JOB_ID']
            except KeyError:
                # If SLURM_JOB_ID is not found, generate a random number as job ID
                self.job_id = str(random.randint(100000, 999999))

        else:
            self.job_id = sim_params.job_id

        sim_params.date_created = pd.Timestamp.now()
        sim_params.species_out_num = grab_species_num(self.run_path, sim_params.species_out)

        # Define input and output paths     
        self.inputs_job = self.data_path + "inputs/" + self.job_id + "/"
        self.inputs_glob = self.data_path + "inputs/glob/"
        self.outputs_job = self.data_path + "outputs/" + self.job_id + "/"

        if os.path.exists(self.inputs_job):
            shutil.rmtree(self.inputs_job)

        if os.path.exists(self.outputs_job):
            shutil.rmtree(self.outputs_job)

        os.makedirs(self.inputs_job)
        os.makedirs(self.outputs_job)

        self.walltimes = {
            "traj_gen_wall": 0,
            "gen_met_wall": 0,
            "bg_chem_wall": 0,
            "ac_perf_wall": 0,
            "emissions_wall": 0,
            "sim_plumes_wall": 0,
            "plume_to_grid_wall": 0,
            "run_cc_wall": 0,
            "run_boxm_wall": 0,
            "gen_outputs_wall": 0,
        }

        self.proctimes = {
            "traj_gen_proc": 0,
            "gen_met_proc": 0,
            "bg_chem_proc": 0,
            "ac_perf_proc": 0,
            "emissions_proc": 0,
            "sim_plumes_proc": 0,
            "plume_to_grid_proc": 0,
            "run_cc_proc": 0,
            "run_boxm_proc": 0,
            "gen_outputs_proc": 0,
        }

        all_params = {
            "fl_params": fl_params,
            "plume_params": plume_params,
            "sim_params": sim_params,
        }

        # Set the model parameters
        self.fl_params = fl_params
        self.plume_params = plume_params
        self.sim_params = sim_params
        self.all_params = all_params

    def eval(self):
        """Run the GPAT model."""

        # Generate formation flight trajectory points
        start_wall_time = time.time()
        start_process_time = time.process_time()
        if self.fl_params.n_ac > 0:
            self.fl = self.traj_gen()
        self.walltimes["traj_gen_wall"] = time.time() - start_wall_time
        self.proctimes["traj_gen_proc"] = time.process_time() - start_process_time

        # Generate meteorological data
        start_wall_time = time.time()
        start_process_time = time.process_time()
        self.met = self.gen_met()
        self.walltimes["gen_met_wall"] = time.time() - start_wall_time
        self.proctimes["gen_met_proc"] = time.process_time() - start_process_time

        # Generate background chemistry data
        start_wall_time = time.time()
        start_process_time = time.process_time()
        self.bg_chem = self.gen_bg_chem()
        self.walltimes["bg_chem_wall"] = time.time() - start_wall_time
        self.proctimes["bg_chem_proc"] = time.process_time() - start_process_time

        # Calculate aircraft performance using PS Model
        start_wall_time = time.time()
        start_process_time = time.process_time()
        if self.fl_params.n_ac > 0:
            self.fl = self.ac_perf()
        self.walltimes["ac_perf_wall"] = time.time() - start_wall_time
        self.proctimes["ac_perf_proc"] = time.process_time() - start_process_time

        # Estimate emissions using Pycontrails Emissions Model
        start_wall_time = time.time()
        start_process_time = time.process_time()
        if self.fl_params.n_ac > 0:
            self.fl = self.emissions()
        self.walltimes["emissions_wall"] = time.time() - start_wall_time
        self.proctimes["emissions_proc"] = time.process_time() - start_process_time

        # Simulate plume dispersion/advection using Pycontrails Dry Advection Model
        start_wall_time = time.time()
        start_process_time = time.process_time()
        if self.fl_params.n_ac > 0:
            self.fl, self.pl = self.sim_plumes()
        self.walltimes["sim_plumes_wall"] = time.time() - start_wall_time
        self.proctimes["sim_plumes_proc"] = time.process_time() - start_process_time

        # Aggregate plumes to an Eulerian grid for photochemical and microphysical processing
        start_wall_time = time.time()
        start_process_time = time.process_time()
        self.emi = self.plume_to_grid()
        self.walltimes["plume_to_grid_wall"] = time.time() - start_wall_time
        self.proctimes["plume_to_grid_proc"] = time.process_time() - start_process_time

        # Run COCIP
        start_wall_time = time.time()
        start_process_time = time.process_time()
        # self.contrail = self.run_cc()
        self.walltimes["run_cc_wall"] = time.time() - start_wall_time
        self.proctimes["run_cc_proc"] = time.process_time() - start_process_time

        # Run BOXM
        start_wall_time = time.time()
        start_process_time = time.process_time()
        self.chem = self.run_boxm()
        self.walltimes["run_boxm_wall"] = time.time() - start_wall_time
        self.proctimes["run_boxm_proc"] = time.process_time() - start_process_time

        # Generate outputs
        start_wall_time = time.time()
        start_process_time = time.process_time()
        self.gen_outputs()
        self.walltimes["gen_outputs_wall"] = time.time() - start_wall_time
        self.proctimes["gen_outputs_proc"] = time.process_time() - start_process_time

    # Model methods
    def traj_gen(self) -> list[Flight]:
        """Generate formation flight trajectory points."""
        fl_params = self.fl_params
        fl = []

        lat0, lon0, alt0 = fl_params.fl0_coords0
        heading = fl_params.fl0_heading
        dist = fl_params.fl0_speed * fl_params.rt_fl.total_seconds()

        # calculate the final coordinates
        geod = Geod(ellps="WGS84")
        lon1, lat1, _ = geod.fwd(lon0, lat0, heading, dist)

        # create flight object for leader flight and resample points according to ts_fl
        df = pd.DataFrame()
        df["longitude"] = [lon0, lon1]
        df["latitude"] = [lat0, lat1]
        df["altitude"] = [alt0, alt0]
        df["time"] = [fl_params.t0_fl, (fl_params.t0_fl + fl_params.rt_fl)]

        ts_fl_min = int(fl_params.ts_fl.total_seconds() / 60)

        fl0 = Flight(df).resample_and_fill(freq=f"{ts_fl_min}min")
        fl0.attrs = {"flight_id": 0, "aircraft_type": fl_params.ac_type}
        mask = (
                    (fl0["latitude"] > self.sim_params.lat_bounds[0] + 0.01) & 
                    (fl0["latitude"] < self.sim_params.lat_bounds[1] - 0.01) &
                    (fl0["longitude"] > self.sim_params.lon_bounds[0] + 0.01) & 
                    (fl0["longitude"] < self.sim_params.lon_bounds[1] - 0.01) &
                    (fl0["altitude"] > self.sim_params.alt_bounds[0]) & 
                    (fl0["altitude"] < self.sim_params.alt_bounds[1])

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
                    (fli["latitude"] > self.sim_params.lat_bounds[0] + 0.01) & 
                    (fli["latitude"] < self.sim_params.lat_bounds[1] - 0.01) &
                    (fli["longitude"] > self.sim_params.lon_bounds[0] + 0.01) & 
                    (fli["longitude"] < self.sim_params.lon_bounds[1] - 0.01) &
                    (fli["altitude"] > self.sim_params.alt_bounds[0]) & 
                    (fli["altitude"] < self.sim_params.alt_bounds[1])
                )
                fli = fli.filter(mask)
                fl.append(fli)

                # Update starting coordinates for next flight
                lon0, lat0, alt0 = lon_dx_dy, lat_dx_dy, alt_dx_dy

        return fl

    def gen_met(self) -> MetDataset:
        """Generate meteorological data."""
        sim_params = self.sim_params

        met = xr.Dataset(
            data_vars={
                "eastward_wind": (("latitude", "longitude", "level", "time"), 
                                  da.full((len(self.lats), len(self.lons), len(self.alts), 
                                           len(self.times)), sim_params.eastward_wind)),
                "northward_wind": (("latitude", "longitude", "level", "time"), 
                                   da.full((len(self.lats), len(self.lons), len(self.alts), 
                                            len(self.times)), sim_params.northward_wind)),
                "lagrangian_tendency_of_air_pressure": (("latitude", "longitude", "level", "time"), 
                    da.full((len(self.lats), len(self.lons), 
                             len(self.alts), len(self.times)), 
                             sim_params.lagrangian_tendency_of_air_pressure)),
                "air_temperature": (("latitude", "longitude", "level", "time"), 
                                    da.zeros((len(self.lats), len(self.lons), len(self.alts), 
                                              len(self.times)))),
            },

            coords={
                "longitude": self.lons, "latitude": self.lats, 
                "level": units.m_to_pl(self.alts), "time": self.times
            }
        )

        met = MetDataset(met)

        month = self.times[0].month

        air_temperature = xr.open_dataarray(
            self.inputs_glob + "air_temperature.nc", engine='netcdf4'
        ).sel(month=month - 1).interp(
            longitude=self.lons, latitude=self.lats, level=self.levels,
            method="linear").broadcast_like(met.data)

        h2o_concs = xr.open_dataarray(
            self.inputs_glob + "h2o_concs.nc", engine='netcdf4'
        ).sel(month=month - 1).interp(
            longitude=self.lons, latitude=self.lats, level=self.levels,
            method="linear").broadcast_like(met.data)

        met.data["air_temperature"] = air_temperature.transpose(
            "latitude", "longitude", "level", "time"
            )

        met.data["H2O"] = h2o_concs.transpose("latitude", "longitude", "level", "time")

        rho_d = met["air_pressure"].data / (constants.R_d * met["air_temperature"].data)

        N_A = 6.022e23  # Avogadro's number

        met.data["specific_humidity"] = met.data["H2O"] * constants.M_d / (N_A * rho_d * 1e-6)

        met.data["relative_humidity"] = thermo.rhi(
            met.data["specific_humidity"],
            met.data["air_temperature"],
            met.data["air_pressure"]
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
                met["latitude"].data.values, 
                met["longitude"].data.values, 
                met["time"].data.values
            ),
        )

        return met

    def gen_bg_chem(self) -> xr.Dataset:
        """Generate background chemistry data."""
        month = self.times[0].month

        bg_chem = xr.open_dataset(
            self.inputs_glob + "species.nc", engine='netcdf4'
        ).sel(month=month - 1)

        for s in [1, 2, 3, 5, 7, 9, 10, 13, 15, 16, 17, 18, 19, 20, 22, 24, 26, 27, 29, 31, 33, 
                  35, 36, 37, 38, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 
                  58, 60, 62, 63, 65, 66, 68, 69, 70, 72, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 
                  85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102, 104, 105, 
                  106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 
                  122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 
                  138, 139, 140, 141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 
                  155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 
                  171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 
                  187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 200, 201, 203, 204, 
                  205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219]:
            bg_chem.bg_chem[:,:,:,s-1] = 0

        bg_chem = bg_chem * 1e09  # convert mixing ratio to ppb
        
        # downselect and interpolate bg_chem to the simulation grid
        return bg_chem.interp(
            longitude=self.lons, latitude=self.lats, level=self.levels
        )

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
        plume_params = self.plume_params
        #met = self.met
        fl = self.fl

        emi_model = Emissions()

        for i, _fli in enumerate(fl):

            # get emissions data
            fl[i] = emi_model.eval(fl[i])

            # Iterate over the columns in the DataFrame
            for column in fl[i].dataframe.columns:
                # Replace NaN values in the column with the value from the previous row
                fl[i].dataframe[column] = fl[i].dataframe[column].fillna(method='ffill')

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
                fl[i][species] = (
                    ei * fl[i]["fuel_flow"] * plume_params.dt_integration.seconds
                )

        return fl

    def sim_plumes(self) -> list[pd.DataFrame]:
        """Simulate plume dispersion/advection using Pycontrails Dry Advection Model."""
        plume_params = self.plume_params
        met = self.met
        fl = self.fl

        # converting from dataclass to dict to pass to DryAdvection model
        plume_params_dict = asdict(plume_params)

        # Create a new dictionary excluding hres_pl and vres_pl to input to dry advection model
        filtered_plume_params_dict = {key: value for key, value in plume_params_dict.items() if 
                                      key not in {"hres_pl", "vres_pl", "n_slices"}}

        dry_adv = DryAdvection(met, **filtered_plume_params_dict)

        pl = []

        for i, fli in enumerate(fl):

            pli = dry_adv.eval(fli)
            pl.append(pli)

            # convert both flights and plumes to dataframes
            fl[i] = fl[i].dataframe
            for column in fl[i].columns:
                # Replace NaN values in the column with the value from the previous row
                fl[i][column] = fl[i][column].fillna(method='ffill')

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
                    "heading",
                    "sigma_yy",
                    "sigma_zz",
                ]
            ],
            on=["flight_id", "waypoint"],
        ).sort_values(by=["time", "flight_id", "waypoint"])

        pl["sin_a"] = np.sin(np.radians(pl["heading"]))
        pl["cos_a"] = np.cos(np.radians(pl["heading"]))
        pl["altitude"] = units.pl_to_m(pl["level"])
        pl["time"] = pl["time"] - self.plume_params.dt_integration

        return fl, pl

    def plume_to_grid(self) -> MetDataset:
        """Aggregate plumes to an Eulerian grid for chemical and physical processing."""
        
        # loop over time and plume property
        emi = xr.DataArray(
            np.zeros((len(self.lons_pl), len(self.lats_pl), len(self.alts), 
                      len(self.times), 9)),
        dims=["longitude", "latitude", "level", "time", "emi_species"],
        coords={
                "longitude": self.lons_pl,
                "latitude": self.lats_pl,
                "level": units.m_to_pl(self.alts),
                "time":  self.times,
                "emi_species": ["NO", "NO2", "CO", "HCHO", "CH3CHO", "C2H4", "C3H6", "C2H2", 
                                "BENZENE"],
            }
        )

        if self.fl_params.n_ac >= 1:

            plume_params = self.plume_params
            #fl = self.fl
            pl = self.pl
            
            # define molar masses of species g/mol
            mm = [30.01, 46.01, 28.01, 30.03, 44.05, 28.05, 42.08, 26.04, 78.11]  # g/mol
            NA = 6.022e23  # Avogadro's number
            bbox = (
                self.lons[0],
                self.lats[0],
                self.lons[-1],
                self.lats[-1],
                self.alts[0],
                self.alts[-1],
            )

            max_age = self.plume_params.max_age  # Define the maximum age for plume waypoints

            # Placeholder for background emissions
            bg_property_data = {property: 0 for property in emi["emi_species"].values}  
            

            
            for t, time in enumerate(pl["time"].unique()):
                print("Processing time: ", time)
                # create geovectordataset to store instantaneous plume data
                plume_time_data = GeoVectorDataset(data=pl.loc[pl["time"] == time])
                calc_continuous(plume_time_data)

                # Find max age plume waypoints
                max_age_segments = plume_time_data.dataframe.loc[
                    plume_time_data.dataframe["age"] == max_age
                    ]

                for p, property in enumerate(emi["emi_species"].values):
                    if property in self.sim_params.species_in:

                        # call contrails_to_hi_res_grid
                        plume_property_data = plume_to_grid(
                            time=time,
                            plumes_t=plume_time_data,
                            var_name=property,
                            spatial_bbox=bbox,
                            spatial_grid_res=plume_params.hres_pl,
                            n_slices=plume_params.n_slices,
                        )

                        # add background emissions mass
                        plume_property_data += bg_property_data[property]

                        # now that aggregation is done, can add max age plumes to 
                        # bg_property_data
                        bg_property_data[property] += max_age_segments[property].sum() / \
                        (len(self.lons) * len(self.lats))

                        # convert mass to density [kg/m^3]
                        density = plume_property_data / (plume_params.vres_pl \
                        * units.latitude_distance_to_m(plume_params.hres_pl) \
                        * units.longitude_distance_to_m(plume_params.hres_pl, 
                                                        (self.lats[0] + self.lats[-1]) / 2))

                        plume = (density / 1E+03) * NA / mm[p] 
                        # [kg/m^3] to [molecules/cm^3]
                        # kg -> g (* 1E+03)
                        # m^3 -> cm^3 (/ 1E+06)

                        # find altitude index for flight level
                        alt = units.m_to_pl(self.fl_params.fl0_coords0[2])

                        # Ensure consistent dimensions
                        plume = plume.reindex(latitude=emi.latitude, longitude=emi.longitude, 
                                              method="nearest")

                        # find time index for time slice covering dt_integration
                        if time < pl["time"].unique()[-1]:
                            next_time = pl["time"].unique()[t + 1]
                            time_slice = slice(time, next_time)
                            emi.loc[:, :, alt, time_slice, property] = plume

                        if time == pl["time"].unique()[-1]:
                            end_time = emi["time"].max()
                            time_slice = slice(time, end_time)
                            emi.loc[:, :, alt, time_slice, property] = plume

        return MetDataset(xr.Dataset({"emi": emi}))

    def run_cc(self) -> xr.Dataset:
        """Run Contrail Model."""
        # met = self.met
        # emi = self.emi

        # return contrail
        pass 
    
    def run_boxm(self) -> xr.Dataset:
        """Run BOXM."""
        # Initialize the box model dataset
        self.init_boxm_ds()
        # Stack the box model dataset
        self.stack()
        # Convert the datasets to netCDF
        self.to_netcdf()
        # Run the box model
        self.do_boxm()
        # Unstack the box model dataset
        self.unstack()

        return self.boxm_ds_unstacked

    def gen_outputs(self):
        """Generate outputs."""

        print("Generating outputs...")
        # Add job runtime for all methods
        self.all_params["walltimes"] = self.walltimes
        self.all_params["proctimes"] = self.proctimes
        
        # Save to pickle file
        with open(f"{self.outputs_job}params_{self.job_id}.pkl", 'wb') as pkl_file:
            pickle.dump(self.all_params, pkl_file)

        # Save fl dataset to pickle file
        if self.fl_params.n_ac > 0:
            self.fl.to_pickle(f"{self.outputs_job}fl_{self.job_id}.pkl")

        # Save pl dataset to pickle file
        if self.fl_params.n_ac > 0:
            self.pl.to_pickle(f"{self.outputs_job}pl_{self.job_id}.pkl")

        # Save the box model dataset to netCDF file
        print("Saving chem dataset to netCDF file...")
        self.chem.to_netcdf(f"{self.outputs_job}chem_{self.job_id}.nc")
        print("Done!")

# Methods for running the box model

    def init_boxm_ds(self):
        """Initialize the box model dataset."""

        self.boxm_ds = xr.merge([self.met.data, self.bg_chem, self.emi.data])

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

        self.boxm_ds = self.boxm_ds.assign_attrs(
            dts=self.sim_params.ts_sim.total_seconds(),
            species_out=self.sim_params.species_out,
            species_out_num = self.sim_params.species_out_num
            )

        self.boxm_ds["J"] = (["time", "level", "longitude", "latitude", "photol_params"], 
                             da.zeros((self.boxm_ds.sizes["time"], 
                                       self.boxm_ds.sizes["level"], 
                                       self.boxm_ds.sizes["longitude"], 
                                       self.boxm_ds.sizes["latitude"], 5)))

        self.boxm_ds["DJ"] = (["time", "level", "longitude", "latitude", "photol_coeffs"], 
                              da.zeros((self.boxm_ds.sizes["time"], 
                                        self.boxm_ds.sizes["level"], 
                                        self.boxm_ds.sizes["longitude"], 
                                        self.boxm_ds.sizes["latitude"], 5)))

        self.boxm_ds["RC"] = (["time", "level", "longitude", "latitude", "therm_coeffs"], 
                              da.zeros((self.boxm_ds.sizes["time"], 
                                        self.boxm_ds.sizes["level"], 
                                        self.boxm_ds.sizes["longitude"], 
                                        self.boxm_ds.sizes["latitude"], 5)))

        self.boxm_ds["Y"] = (["time", "level", "longitude", "latitude", "species_out"], 
                             da.zeros((self.boxm_ds.sizes["time"], 
                                       self.boxm_ds.sizes["level"], 
                                       self.boxm_ds.sizes["longitude"], 
                                       self.boxm_ds.sizes["latitude"], 
                                       len(self.sim_params.species_out))))

    def stack(self):
        """Stack boxm_ds to flatten and get cell numbers out."""

        # stack datasets to get cell index for fortran
        self.boxm_ds_stacked = self.boxm_ds.stack(
            {"cell": ["level", "longitude", "latitude"]}
        )

        self.boxm_ds_stacked = self.boxm_ds_stacked.reset_index("cell")

    def to_netcdf(self):
        """Convert the met, bg_chem, and emi datasets to boxm_ds.nc for use in the box model."""

        # Delete any existing netCDF files
        if pathlib.Path(f"{self.inputs_job}/boxm_ds.nc").exists():
            print("deleting boxm_ds.nc")
            pathlib.Path(f"{self.inputs_job}/boxm_ds.nc").unlink()

        # Convert DataFrames to Datasets and write to netCDF
        self.boxm_ds_stacked.to_netcdf(f"{self.inputs_job}/boxm_ds.nc", mode="w")

    def do_boxm(self):
        """Run the box model in fortran using subprocess."""

        # Run the box model
        subprocess.call(
            [self.run_path + "boxm", self.job_id], 
        )

        # open nc file
        self.boxm_ds = xr.open_dataset(f"{self.inputs_job}/boxm_ds.nc")

    def unstack(self):
        """Unstack the box model dataset."""
        print("Chunking the dataset")
        # Convert the dataset to a Dask dataset
        self.boxm_ds = self.boxm_ds.chunk({'cell': 100})  # Adjust chunk size based on your mem

        print("Set coords")
        # Convert 'level', 'lat', and 'lon' to coordinates
        self.boxm_ds_unstacked = self.boxm_ds.set_coords(['level', 'longitude', 'latitude'])

        print("Set index")
        # Create a multi-index for the 'cell' dimension
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.set_index(
            cell=['level', 'longitude', 'latitude']
            )

        print("Unstack the dataset")
        # Unstack the dataset
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.unstack("cell")

        print("Compute the result")
        # # Compute the result to trigger the lazy evaluation
        #self.boxm_ds_unstacked = self.boxm_ds_unstacked.compute()


# Functions used in GPAT Model
def grab_species_num(run_path, species_out: np.array) -> np.array:
    """Grab the species numbers for the species of interest in output."""
    # Read species names from the file into a list
    with open(f'{run_path}species_num.txt') as file:
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
        np.concatenate([[heading[0]], heading]) if len(heading) > 0 else [np.nan], 
        index=group.index
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
            
            previous_time = pl_df["time"].unique()[ts-1]
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
                grid_concs_over_zero = (grid_concs > 0)
                grid_concs_over_zero = grid_concs.where(grid_concs_over_zero, drop=True)
                               
                grid_mass = grid_concs_over_zero \
                    * chem_ds["M"].sel(time=t) \
                    * 1e-9 \
                    * (mm[s] / NA) \
                    * params.loc["vres_sim"] \
                    * units.latitude_distance_to_m(params.loc["hres_sim"]) \
                    * units.longitude_distance_to_m(
                        params.loc["hres_sim"], 
                        (params.loc["lat_bounds"][0] + params.loc["lat_bounds"][1]) / 2
                        ) \
                    * 1E+03  # convert to kg/m^3

                total_grid_mass = grid_mass.sum().item()

                percent_mass_conserved = total_grid_mass / total_vector_mass * 100
                
            # Append the percentage to the list in the dictionary
            vecmass[emi_species].append(total_vector_mass)
            gridmass[emi_species].append(total_grid_mass)
            mc[emi_species].append(percent_mass_conserved)

    # convert the dictionary to a DataFrame
    vecmass = pd.DataFrame(vecmass, index=pl_df["time"].unique()[:-1], 
                           columns=chem_ds["emi_species"].values.tolist())
    gridmass = pd.DataFrame(gridmass, index=pl_df["time"].unique()[:-1], 
                            columns=chem_ds["emi_species"].values.tolist())
    mc = pd.DataFrame(mc, index=pl_df["time"].unique()[:-1], 
                      columns=chem_ds["emi_species"].values.tolist())

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
        runtime = int((end_time - start_time) / np.timedelta64(1, 'D')) % 365
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
    for s in ["NO2", "NO", "O3", "CO", "CH4", "HCHO", "CH3CHO", "CH3COCH3",
                        "C2H6", "C2H4", "C3H8", "C3H6", "C2H2", "NC4H10", "TBUT2ENE",
                        "BENZENE", "TOLUENE", "OXYL", "C5H8", "H2O2", "HNO3", "C2H5CHO",
                        "CH3OH", "MEK", "CH3OOH", "PAN", "MPAN"]:
        
        boxm_input.write(f"{cell_chem_ds["bg_chem"].sel(species=s).item()}\n")
        
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

    sza_df = pd.read_csv(f"{data_path}outputs/{job_id}/ZEN.OUT", header=0,
                        names=['TIME', 'ZEN'], dtype=np.float64)
        
    J_df = pd.read_csv(f"{data_path}outputs/{job_id}/J.OUT", header=0,
                        names=['TIME', 'J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 
                               'J10', 'J11', 'J12', 'J13','J14', 'J15', 'J16', 'J17', 'J18', 
                               'J19', 'J20', 'J21', 'J22', 'J23', 'J24', 'J25', 'J26', 'J27', 
                               'J28', 'J29', 'J30', 'J31', 'J32', 'J33', 'J34', 'J35', 'J36', 
                               'J37', 'J38', 'J39', 'J40', 'J41', 'J42', 'J43', 'J44', 'J45', 
                               'J46', 'J47', 'J48', 'J49', 'J50'], dtype=np.float64)

    DJ_df = pd.read_csv(f"{data_path}outputs/{job_id}/DJ.OUT", header=0,
                            names=['TIME', 'DJ1', 'DJ2', 'DJ3', 'DJ4', 'DJ5', 'DJ6', 'DJ7', 'DJ8', 
                                   'DJ9', 'DJ10', 'DJ11', 'DJ12', 'DJ13','DJ14', 'DJ15', 'DJ16', 
                                   'DJ17', 'DJ18', 'DJ19', 'DJ20', 'DJ21', 'DJ22', 'DJ23', 'DJ24', 
                                   'DJ25', 'DJ26', 'DJ27', 'DJ28', 'DJ29', 'DJ30', 'DJ31', 'DJ32', 
                                   'DJ33', 'DJ34', 'DJ35', 'DJ36', 'DJ37', 'DJ38', 'DJ39', 'DJ40', 
                                   'DJ41', 'DJ42', 'DJ43', 'DJ44', 'DJ45', 'DJ46', 'DJ47', 'DJ48', 
                                   'DJ49', 'DJ50'], dtype=np.float64)

    RC_df = pd.read_csv(f"{data_path}outputs/{job_id}/RC.OUT", header=0,
                        names=['TIME', 'RC1', 'RC2', 'RC3', 'RC4', 'RC5', 'RC6', 'RC7', 'RC8', 
                               'RC9', 'RC10', 'RC11', 'RC12', 'RC13','RC14', 'RC15', 'RC16', 
                               'RC17', 'RC18', 'RC19', 'RC20', 'RC21', 'RC22', 'RC23', 'RC24', 
                               'RC25', 'RC26', 'RC27', 'RC28', 'RC29', 'RC30', 'RC31', 'RC32', 
                               'RC33', 'RC34', 'RC35', 'RC36', 'RC37', 'RC38', 'RC39', 'RC40', 
                               'RC41', 'RC42', 'RC43', 'RC44', 'RC45', 'RC46', 'RC47', 'RC48', 
                               'RC49', 'RC50'], dtype=np.float64)

    # get species names
    header_names = ["TIME", *list(cell_chem_ds["species"].values)]

    Y_df = pd.read_csv(f"{data_path}/outputs/{job_id}/Y.OUT", header=0,
                            names=header_names, dtype=np.float64) 
    
    # # Update the chem_ds_stacked with the new data
    # Update zen data
    cell_chem_ds["sza_orig"] = (["time"], da.zeros(cell_chem_ds.sizes["time"]))
    cell_chem_ds["sza_orig"].loc[:] = sza_df["ZEN"].values * np.pi / 180

    cell_chem_ds["J_orig"] = (["time", "photol_params"], 
                              da.zeros((cell_chem_ds.sizes["time"], 5)))
    for pp, photol_params in enumerate(J_df.columns[1:6]):
        cell_chem_ds["J_orig"].loc[:, pp] = J_df[photol_params].values

    cell_chem_ds["DJ_orig"] = (["time", "photol_coeffs"], 
                               da.zeros((cell_chem_ds.sizes["time"], 5)))
    for pc, photol_coeffs in enumerate(DJ_df.columns[1:6]):
        cell_chem_ds["DJ_orig"].loc[:, pc] = DJ_df[photol_coeffs].values

    cell_chem_ds["RC_orig"] = (["time", "therm_coeffs"], 
                               da.zeros((cell_chem_ds.sizes["time"], 5)))
    for tc, therm_coeffs in enumerate(RC_df.columns[1:6]):
        cell_chem_ds["RC_orig"].loc[:, tc] = RC_df[therm_coeffs].values
        
    cell_chem_ds["Y_orig"] = (["time", "species_out"], 
                              da.zeros((cell_chem_ds.sizes["time"], 
                                        cell_chem_ds.sizes["species_out"])))
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
    parser.add_argument("--shear", type=float, help="Wind shear in 1/s")
    parser.add_argument("--hres_pl", type=float, 
                        help="Horizontal resolution of the plume in degrees")
    parser.add_argument("--vres_pl", type=float, 
                        help="Vertical resolution of the plume in meters")
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
    parser.add_argument("--lagrangian_tendency_of_air_pressure", type=float, 
                        help="Lagrangian tendency of air pressure in m/s")
    parser.add_argument("--species_in", type=str, help="Input species (comma-separated)")
    parser.add_argument("--species_out", type=str, help="Output species (comma-separated)")
    parser.add_argument("--gpat_path", type=str, help="GPAT directory")
    parser.add_argument("--job_id", type=str, help="Job ID")
    parser.add_argument("--run_gpat", action='store_true', help="Run the GPAT model")
       
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
        params.fl0_coords0 = tuple(map(float, args.fl0_coords0.split(',')))
    if args.sep_dist:
        params.sep_dist = tuple(map(float, args.sep_dist.split(',')))
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
    if args.shear:
        params.shear = args.shear
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
        params.lat_bounds = tuple(map(float, args.lat_bounds.split(',')))
    if args.lon_bounds:
        params.lon_bounds = tuple(map(float, args.lon_bounds.split(',')))
    if args.alt_bounds:
        params.alt_bounds = tuple(map(float, args.alt_bounds.split(',')))
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
        params.species_in = tuple(args.species_in.split(','))
    if args.species_out:
        params.species_out = tuple(args.species_out.split(','))
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