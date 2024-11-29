"""Gridded Plume Analysis Tool (GPAT).

Simulate aircraft trajectories, estimate aircraft performance, fuel burn and emissions.

Plot associated aircraft exhaust plumes, subject to Gaussian dispersion and advection. Aggregate plumes to an Eulerian grid for photochemical and microphysical processing."""

import os
import random
import numpy as np
import pandas as pd
import xarray as xr
import dask.array as da
import yaml
import pickle
import time
import shutil
from pyproj import Geod
import scipy.stats as stats
import matplotlib.pyplot as plt
import argparse
from matplotlib.animation import FuncAnimation, PillowWriter
import subprocess
from pycontrails.core import Flight, GeoVectorDataset, MetDataArray, MetDataset, models
from pycontrails.core.models import Model, ModelParams
from pycontrails.models.emissions import Emissions
from pycontrails.models.ps_model import PSFlight
from pycontrails.models.dry_advection import DryAdvection
from pycontrails.models.gpat.plume_to_grid import plume_to_grid
from pycontrails.models.cocip import contrails_to_hi_res_grid
from pycontrails.physics import geo, thermo, units, constants
from dataclasses import dataclass, asdict, fields, is_dataclass
from distutils.util import strtobool
from typing import Tuple
import pathlib

@dataclass
class SimParams():
    """Default simulation parameters"""
    t0_sim: pd.Timestamp = pd.to_datetime("2022-01-20 12:00:00") # simulation start time
    rt_sim: pd.Timedelta = pd.Timedelta(hours=120) # simulation runtime
    ts_sim: pd.Timedelta = pd.Timedelta(seconds=20) # simulation time step
    lat_bounds: Tuple[float, float] = (0.0, 1.0) # lat bounds [deg]
    lon_bounds: Tuple[float, float] = (0.0, 1.0) # lon bounds [deg]
    alt_bounds: Tuple[float, float] = (12000, 13000) # alt bounds [m]
    hres_sim: float = 0.01 # horizontal resolution [deg]
    vres_sim: float = 500 # vertical resolution [m]
    eastward_wind: float = 0.0 # m/s
    northward_wind: float = 0.0 # m/s
    lagrangian_tendency_of_air_pressure: float = 0.0 # m/s
    species_in: Tuple = ("NO")
    species_out: Tuple = ("O3", "NO2", "NO",
                        "NO3", "N2O5", "HNO3",
                        "HONO", "HO2", "OH",
                        "H2O2", "H2O", "CO",
                        "CH4", "C2H6", "C3H8",
                        "C2H4", "C3H6")
    gpat_path: str = None   
    job_id: str = None
    run_gpat: bool = False
    date_created: pd.Timestamp = None
    species_out_num: Tuple = None

class GPAT_clean(Model):
    """Gridded Plume Analysis Tool (GPAT).

    Simulate photochemical and microphysical processing of atmosphere.

    Parameters
    ----------
    sim_params : SimParams
        Simulation parameters.
    """

    name = "GPAT_clean"
    long_name = "Gridded Plume Analysis Tool no emission version"
    #default_params = (FlParams, PlumeParams, SimParams)

    def __init__(
            self,
            sim_params: SimParams
            ):
        super().__init__()

        # Generate the grid
        self.lats = np.arange(
            sim_params.lat_bounds[0], sim_params.lat_bounds[1] + sim_params.hres_sim, sim_params.hres_sim
        )
        self.lons = np.arange(
            sim_params.lon_bounds[0], sim_params.lon_bounds[1] + sim_params.hres_sim, sim_params.hres_sim
        )
        self.alts = np.arange(
            sim_params.alt_bounds[0], sim_params.alt_bounds[1] + sim_params.vres_sim, sim_params.vres_sim
        )
        self.levels = units.m_to_pl(self.alts)

        self.times = pd.date_range(
            start=sim_params.t0_sim,
            end=sim_params.t0_sim + sim_params.rt_sim,
            freq=sim_params.ts_sim,
        )

        self.total_volume = (
            (self.lats[-1] - self.lats[0])
            * (self.lons[-1] - self.lons[0])
            * (self.alts[-1] - self.alts[0])
        )

        if sim_params.gpat_path is None:
            self.path = os.environ['PYCONTRAILSDIR'] + "models/gpat/"

        else:
            self.path = sim_params.gpat_path

        if sim_params.job_id is None:
            try:
                self.job_id = os.environ['SLURM_JOB_ID']
            except KeyError:
                # If SLURM_JOB_ID is not found, generate a random number as job ID
                self.job_id = str(random.randint(100000, 999999))

        else:
            self.job_id = sim_params.job_id

        sim_params.date_created = pd.Timestamp.now()
        sim_params.species_out_num = grab_species_num(sim_params.species_out)

        # Define input and output paths     
        self.inputs_job = self.path + "inputs/" + self.job_id + "/"
        self.inputs_glob = self.path + "inputs/glob/"
        self.outputs_job = self.path + "outputs/" + self.job_id + "/"

        if os.path.exists(self.inputs_job):
            shutil.rmtree(self.inputs_job)

        if os.path.exists(self.outputs_job):
            shutil.rmtree(self.outputs_job)

        os.makedirs(self.inputs_job)
        os.makedirs(self.outputs_job)

        self.walltimes = {
            "gen_met_wall": 0,
            "bg_chem_wall": 0,
            "run_cc_wall": 0,
            "run_boxm_wall": 0,
            "gen_outputs_wall": 0,
        }

        self.proctimes = {
            "gen_met_proc": 0,
            "bg_chem_proc": 0,
            "run_cc_proc": 0,
            "run_boxm_proc": 0,
            "gen_outputs_proc": 0,
        }

        all_params = {
            "sim_params": sim_params,
        }

        # Set the model parameters
        self.sim_params = sim_params
        self.all_params = all_params

    def eval(self):
        """Run the GPAT model."""
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
    def gen_met(self) -> MetDataset:
        """Generate meteorological data."""
        sim_params = self.sim_params

        met = xr.Dataset(
            data_vars={
                "eastward_wind": (("latitude", "longitude", "level", "time"), da.full((len(self.lats), len(self.lons), len(self.alts), len(self.times)), sim_params.eastward_wind)),
                "northward_wind": (("latitude", "longitude", "level", "time"), da.full((len(self.lats), len(self.lons), len(self.alts), len(self.times)), sim_params.northward_wind)),
                "lagrangian_tendency_of_air_pressure": (("latitude", "longitude", "level", "time"), da.full((len(self.lats), len(self.lons), len(self.alts), len(self.times)), sim_params.lagrangian_tendency_of_air_pressure)),
                "air_temperature": (("latitude", "longitude", "level", "time"), da.zeros((len(self.lats), len(self.lons), len(self.alts), len(self.times)))),
            },

            coords={
                "longitude": self.lons, "latitude": self.lats, "level": units.m_to_pl(self.alts), "time": self.times
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

        met.data["air_temperature"] = air_temperature.transpose("latitude", "longitude", "level", "time")

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
                met["latitude"].data.values, met["longitude"].data.values, met["time"].data.values
            ),
        )

        return met

    def gen_bg_chem(self) -> xr.Dataset:
        """Generate background chemistry data."""
        month = self.times[0].month

        bg_chem = xr.open_dataset(
            self.inputs_glob + "species.nc", engine='netcdf4'
        ).sel(month=month - 1)
        print(bg_chem["bg_chem"])
        # for s in [1, 2, 3, 5, 7, 9, 10, 13, 15, 16, 17, 18, 19, 20, 22, 24, 26, 27, 29, 31, 33, 35, 36, 37, 38, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 62, 63, 65, 66, 68, 69, 70, 72, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 200, 201, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219]:
            
        #     bg_chem.isel(species=s-1)[:,:,:] = 0
        #     print(bg_chem.isel(species=s-1)[:,:,:])

        for s in [1, 2, 3, 5, 7, 9, 10, 13, 15, 16, 17, 18, 19, 20, 22, 24, 26, 27, 29, 31, 33, 35, 36, 37, 38, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 62, 63, 65, 66, 68, 69, 70, 72, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 200, 201, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219]:
            bg_chem.bg_chem[:,:,:,s-1] = 0

        bg_chem = bg_chem * 1e09  # convert mixing ratio to ppb
        
        # downselect and interpolate bg_chem to the simulation grid
        bg_chem = bg_chem.interp(
            longitude=self.lons, latitude=self.lats, level=self.levels
        )

        return bg_chem

    def run_cc(self) -> xr.Dataset:
        """Run Contrail Model."""
        met = self.met
        emi = self.emi

        return contrail

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

        chem = self.boxm_ds_unstacked

        return chem

    def gen_outputs(self):

        print("Generating outputs...")
        # Add job runtime for all methods
        self.all_params["walltimes"] = self.walltimes
        self.all_params["proctimes"] = self.proctimes
        
        # Save to pickle file
        with open(f"{self.outputs_job}params_{self.job_id}.pkl", 'wb') as pkl_file:
            pickle.dump(self.all_params, pkl_file)

        # # Assign species_out names to coord
        # self.chem = self.chem.assign_coords(species_out=self.sim_params.species_out)

        # Save the box model dataset to netCDF file
        print("Saving chem dataset to netCDF file...")
        self.chem.to_netcdf(f"{self.outputs_job}chem_{self.job_id}.nc")
        print("Done!")

    # Methods for running the box model
    def init_boxm_ds(self):

        self.boxm_ds = xr.merge([self.met.data, self.bg_chem])

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

        self.boxm_ds["J"] = (["time", "level", "longitude", "latitude", "photol_params"], da.zeros((self.boxm_ds.sizes["time"], self.boxm_ds.sizes["level"], self.boxm_ds.sizes["longitude"], self.boxm_ds.sizes["latitude"], 5)))

        self.boxm_ds["DJ"] = (["time", "level", "longitude", "latitude", "photol_coeffs"], da.zeros((self.boxm_ds.sizes["time"], self.boxm_ds.sizes["level"], self.boxm_ds.sizes["longitude"], self.boxm_ds.sizes["latitude"], 5)))

        self.boxm_ds["RC"] = (["time", "level", "longitude", "latitude", "therm_coeffs"], da.zeros((self.boxm_ds.sizes["time"], self.boxm_ds.sizes["level"], self.boxm_ds.sizes["longitude"], self.boxm_ds.sizes["latitude"], 5)))

        self.boxm_ds["Y"] = (["time", "level", "longitude", "latitude", "species_out"], da.zeros((self.boxm_ds.sizes["time"], self.boxm_ds.sizes["level"], self.boxm_ds.sizes["longitude"], self.boxm_ds.sizes["latitude"], len(self.sim_params.species_out))))

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
            [self.path + "boxm_clean", self.job_id], 
        )

        # open nc file
        self.boxm_ds = xr.open_dataset(f"{self.inputs_job}/boxm_ds.nc")

    def unstack(self):
        """Unstack the box model dataset."""
        print("Chunking the dataset")
        # Convert the dataset to a Dask dataset
        self.boxm_ds = self.boxm_ds.chunk({'cell': 100})  # Adjust chunk size based on your memory

        print("Set coords")
        # Convert 'level', 'lat', and 'lon' to coordinates
        self.boxm_ds_unstacked = self.boxm_ds.set_coords(['level', 'longitude', 'latitude'])

        print("Set index")
        # Create a multi-index for the 'cell' dimension
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.set_index(cell=['level', 'longitude', 'latitude'])

        print("Unstack the dataset")
        # Unstack the dataset
        self.boxm_ds_unstacked = self.boxm_ds_unstacked.unstack("cell")

        print("Compute the result")
        # # Compute the result to trigger the lazy evaluation
        #self.boxm_ds_unstacked = self.boxm_ds_unstacked.compute()


# Functions used in GPAT Model
def grab_species_num(species_out: np.array) -> np.array:
    """Grab the species numbers for the species of interest in output."""
    # Read species names from the file into a list
    with open('species_num.txt', 'r') as file:
        species_list = [line.strip() for line in file]

    # Create a dictionary mapping species names to their line numbers
    species_dict = {species: index for index, species in enumerate(species_list)}

    return np.array([species_dict[species] + 1 for species in species_out])

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


### Functions for post-processing
def create_jobs_df(outputs_dir):
    jobs = []

    for job_id in os.listdir(outputs_dir):
        job_dir = os.path.join(outputs_dir, job_id)
        if os.path.isdir(job_dir):
            # Check if the params file exists in the subdirectory
            params_file = os.path.join(job_dir, f"params_{job_id}.pkl")
            if os.path.isfile(params_file):
                params = pd.read_pickle(params_file)
                
                # Flatten the dictionary
                data_dict = {}
                for outer_key, inner_dc in params.items():
                    if is_dataclass(inner_dc):
                        inner_dict = asdict(inner_dc)
                    else:
                        inner_dict = inner_dc
                    for inner_key, inner_value in inner_dict.items():
                        data_dict[f"{inner_key}"] = inner_value

                # Check if n_ac > 0 and verify the existence of fl and pl files
                if data_dict.get("n_ac", 0) > 0:
                    expected_files = [
                        f"params_{job_id}.pkl",
                        f"fl_{job_id}.pkl",
                        f"pl_{job_id}.pkl",
                        f"chem_{job_id}.nc"
                    ]
                else:
                    expected_files = [
                        f"params_{job_id}.pkl",
                        f"chem_{job_id}.nc"
                    ]

                if all(os.path.isfile(os.path.join(job_dir, file)) for file in expected_files):
                    jobs.append(data_dict)

    jobs_df = pd.DataFrame(jobs)
    jobs_df = jobs_df.set_index("job_id")

    return jobs_df

def filter_jobs_df(jobs_df, criteria):
    filtered_df = jobs_df.copy()
    
    for key, value in criteria.items():
        if isinstance(value, tuple) and len(value) == 2:
            # Range filter
            filtered_df = filtered_df[(filtered_df[key] >= value[0]) & (filtered_df[key] <= value[1])]
        elif key == "job_id":
            if isinstance(value, list):
                # Combine the list of strings into a single regex pattern
                pattern = '|'.join(value)
                filtered_df = filtered_df[filtered_df.index.str.contains(pattern)]
            else:
                filtered_df = filtered_df[filtered_df.index.str.contains(value)]
        else:
            # Exact match filter
            filtered_df = filtered_df[filtered_df[key] == value]

    return filtered_df

def load_chem_ds(job_ids, outputs_dir):
    chemistry_data = []
    for job_id in job_ids:
        ds = xr.open_dataset(outputs_dir + job_id + "/chem_" + job_id + ".nc")
        ds = ds.expand_dims(job_id=[job_id])
        chemistry_data.append(ds)
    #return xr.concat(chemistry_data, dim="job_id")
    return chemistry_data

def load_chem_da(job_ids, outputs_dir, property):
    chemistry_data = []
    for job_id in job_ids:
        ds = xr.open_dataset(outputs_dir + job_id + "/chem_" + job_id + ".nc")
        ds = ds.expand_dims(job_id=[job_id])
        chemistry_data.append(ds[property])
    return chemistry_data


# Data visualisation
def plot_heatmap(job_id, jobs_df, fl_df, pl_df, chem_ds, **plot_params):
    fig1, ax1 = plt.subplots()
    ax1.set_xticks(np.arange(chem_ds["longitude"][0], chem_ds["longitude"][-1], 0.05))
    ax1.set_yticks(np.arange(chem_ds["latitude"][0], chem_ds["latitude"][-1], 0.05))

    params = jobs_df.loc[job_id]
    fl_df_job = fl_df.loc[job_id]
    pl_df_job = pl_df.loc[job_id]
    chem_ds_job = chem_ds.sel(job_id=job_id)

    print(f"ts: {plot_params['ts']}")
    print(f"plume time: {pl_df_job['time'].iloc[plot_params['ts']]}")
    print(f"chem time: {chem_ds_job["time"][plot_params["ts"]].item()}")

    print(chem_ds_job)
    # Plot the heatmap
    heatmap_data = (
        chem_ds_job["Y"].sel(species_out=plot_params['property'], time=pl_df_job["time"].iloc[plot_params["ts"]])
        .sel(level=178.6, method="nearest")
        .transpose("latitude", "longitude")
    )
    heatmap_data.plot(ax=ax1, cmap="summer")  # You can choose a colormap of your preference

    scat_fl = ax1.scatter(
        fl_df_job["longitude"].loc[fl_df_job["time"] == pl_df_job["time"].iloc[plot_params["ts"]]],
        fl_df_job["latitude"].loc[fl_df_job["time"] == pl_df_job["time"].iloc[plot_params["ts"]]],
        s=5,
        c="red",
        label="Flight path",
    )

    scat_pl = ax1.scatter(
        pl_df_job["longitude"].loc[pl_df_job["time"] == pl_df_job["time"].iloc[plot_params["ts"]]],
        pl_df_job["latitude"].loc[pl_df_job["time"] == pl_df_job["time"].iloc[plot_params["ts"]]],
        s=10e-2 * pl_df_job["width"].loc[pl_df_job["time"] == pl_df_job["time"].iloc[plot_params["ts"]]],
        c="blue",
        label="Plume evolution",
    )

    ax1.legend(loc="upper left")
    ax1.set_xlim([params["lon_bounds"][0], params["lon_bounds"][1]])
    ax1.set_ylim([params["lat_bounds"][0], params["lat_bounds"][1]])
    plt.grid()
    plt.show()

    # print(boxm_da)

    # times = boxm_da["time"].values
    # times_resampled = pd.to_datetime(times).to_series().resample(resample_freq).asfreq().dropna().index

    # print(f"New number of frames: {len(times_resampled)}")

    # def heatmap_func(t):
    #     ax.cla()
    #     ax.set_title(t)

    #     boxm_da.sel(time=t).transpose("latitude", "longitude").plot(
    #         ax=ax, cbar_kwargs={"cax": cbar_ax}, add_colorbar=True, vmin=boxm_da.min(), vmax=boxm_da.max()
    #     )

    # anim = FuncAnimation(fig, heatmap_func, frames=times_resampled, blit=False)

    # filename = pathlib.Path(self.outputs_plots + var1 + "_" + var2 + ".gif")

    # anim.save(filename, dpi=300, writer=PillowWriter(fps=8))

def anim_chem(job_id, jobs_df, fl_df, pl_df, chem_ds, var1, var2, level, resample_freq='4min'):
    """Animate the chemical concentrations with plume vector data."""
    fig, (ax, cbar_ax) = plt.subplots(
        1, 2, gridspec_kw={"width_ratios": (0.9, 0.05), "wspace": 0.2}, figsize=(12, 8)
    )

    params = jobs_df.loc[job_id]
    fl_df_job = fl_df.loc[job_id]
    pl_df_job = pl_df.loc[job_id]
    chem_ds_job = chem_ds.sel(job_id=job_id, time=chem_ds.time[0:1000])

    if var1 == "Y":
        boxm_da = chem_ds_job[var1].sel(species_out=var2).sel(level=level, method="nearest")

    if var1 == "emi":
        boxm_da = chem_ds_job[var1].sel(emi_species=var2).sel(level=level, method="nearest")

    if var1 == "J":
        boxm_da = chem_ds_job[var1].sel(photol_params=var2).sel(level=level, method="nearest")

    if var1 == "DJ":
        boxm_da = chem_ds_job[var1].sel(photol_coeffs=var2).sel(level=level, method="nearest")

    if var1 == "RC":
        boxm_da = chem_ds_job[var1].sel(therm_coeffs=var2).sel(level=level, method="nearest")

    print(boxm_da)

    times = boxm_da["time"].values
    times_resampled = pd.to_datetime(times).to_series().resample(resample_freq).asfreq().dropna().index

    print(f"New number of frames: {len(times_resampled)}")

    # Initialize the first frame to set up the colorbar
    initial_frame = times_resampled[0]
    heatmap_data = boxm_da.sel(time=initial_frame).transpose("latitude", "longitude")
    heatmap = heatmap_data.plot(ax=ax, cmap="Blues", add_colorbar=False, vmin=boxm_da.min(), vmax=boxm_da.max())
    cbar = plt.colorbar(heatmap, cax=cbar_ax)

    def heatmap_func(t):
        ax.cla()
        ax.set_title(t)

        heatmap_data = boxm_da.sel(time=t).transpose("latitude", "longitude")
        heatmap = heatmap_data.plot(ax=ax, cmap="Blues", add_colorbar=False, vmin=boxm_da.min(), vmax=boxm_da.max())#, #)

        # Update the plume vector data
        fl_data = fl_df_job[fl_df_job["time"] == t]
        pl_data = pl_df_job[pl_df_job["time"] == t]

        # Plot the plume vector data
        scat_fl = ax.scatter(
            fl_data["longitude"],
            fl_data["latitude"],
            s=5,
            c="red",
            label="Flight path",
        )

        # scat_pl = ax.scatter(
        #     pl_data["longitude"],
        #     pl_data["latitude"],
        #     s=10e-2 * pl_data["width"],
        #     c="blue",
        #     label="Plume evolution",
        # )

        ax.legend(loc="upper left")
        ax.set_xlim([params["lon_bounds"][0], params["lon_bounds"][1]])
        ax.set_ylim([params["lat_bounds"][0], params["lat_bounds"][1]])

    anim = FuncAnimation(fig, heatmap_func, frames=times_resampled, blit=False)

    filename = pathlib.Path(f"{var1}_{var2}_{job_id}.gif")

    anim.save(filename, dpi=300, writer=PillowWriter(fps=10))

    plt.show()


# Validation
def mc_test(job_id, jobs_df, fl_df, pl_df, chem_ds):
    """Check if mass is conserved in the box model."""

    params = jobs_df.loc[job_id]
    fl_df_job = fl_df.loc[job_id]
    pl_df_job = pl_df.loc[job_id]
    chem_ds_job = chem_ds.sel(job_id=job_id)


    # Initialize the dictionary
    vecmass = {emi_species: [] for emi_species in chem_ds_job["emi_species"].values.tolist()}
    gridmass = {emi_species: [] for emi_species in chem_ds_job["emi_species"].values.tolist()}
    mc = {emi_species: [] for emi_species in chem_ds_job["emi_species"].values.tolist()}

    # Constants
    mm = [30.01, 46.01, 28.01, 30.03, 44.05, 28.05, 42.08, 26.04, 78.11]  # g/mol
    NA = 6.022e23  # Avogadro's number

    for s, emi_species in enumerate(chem_ds_job["emi_species"].values):
        
        max_fl_time = fl_df_job["time"].max()

        for ts, time in enumerate(pl_df_job["time"].unique()[:-1]):
            if ts == 0:
                total_vector_mass = 0
                total_grid_mass = 0
                percent_mass_conserved = 0
                vecmass[emi_species].append(total_vector_mass)
                gridmass[emi_species].append(total_grid_mass)
                mc[emi_species].append(percent_mass_conserved)
                continue
            
            previous_time = pl_df_job["time"].unique()[ts-1]
            fl_snapshot = fl_df_job[fl_df_job["time"] == previous_time]

            if time <= max_fl_time:
                # Accumulate vector mass for all flights
                vector_mass = fl_snapshot[emi_species]

                total_vector_mass += vector_mass.sum()

            # Grab plume mass from grid data
            grid_concs = chem_ds_job["emi"].sel(emi_species=emi_species, time=time).sel(level=178.6, method="nearest")

            if (grid_concs == 0).all():
                pass
            else:
                grid_concs_over_zero = grid_concs.where(grid_concs > 0, drop=True)

                grid_mass = grid_concs_over_zero \
                    * chem_ds_job["M"].sel(time=time).sel(level=178.6, method="nearest") \
                    * 1e-9 \
                    * (mm[s] / NA) \
                    * params.loc["vres_sim"] \
                    * units.latitude_distance_to_m(params.loc["hres_sim"]) \
                    * units.longitude_distance_to_m(params.loc["hres_sim"], (params.loc["lat_bounds"][0] + params.loc["lat_bounds"][1]) / 2) \
                    * 1E+03  # convert to kg/m^3

                total_grid_mass = grid_mass.sum().item()

                percent_mass_conserved = total_grid_mass / total_vector_mass * 100
                
            # Append the percentage to the list in the dictionary
            vecmass[emi_species].append(total_vector_mass)
            gridmass[emi_species].append(total_grid_mass)
            mc[emi_species].append(percent_mass_conserved)

    # convert the dictionary to a DataFrame
    vecmass = pd.DataFrame(vecmass, index=pl_df_job["time"].unique()[:-1], columns=chem_ds_job["emi_species"].values.tolist())
    gridmass = pd.DataFrame(gridmass, index=pl_df_job["time"].unique()[:-1], columns=chem_ds_job["emi_species"].values.tolist())
    mc = pd.DataFrame(mc, index=pl_df_job["time"].unique()[:-1], columns=chem_ds_job["emi_species"].values.tolist())

    # Save the mass conservation data to a pickle file
    pd.to_pickle(mc, f"outputs/{job_id}/mc_{job_id}.pkl")

    return vecmass, gridmass, mc

def boxm_test(path, job_id, cell, chem_ds):
    """Run the box model for selected cells and job_id."""

    chem_ds_stacked = chem_ds.stack(
            {"cell": ["level", "longitude", "latitude"]}
        )
    chem_ds_stacked = chem_ds_stacked.reset_index("cell")

    chem_ds_stacked = chem_ds_stacked.assign_coords(species_out=chem_ds_stacked.attrs["species_out"])

    cell_chem_ds = chem_ds_stacked.sel(job_id=job_id, cell=cell)

    # create input file for original boxm
    gen_boxm_orig_input(cell_chem_ds, job_id)

    gen_zen_file(cell_chem_ds, job_id)

    gen_emi_file(cell_chem_ds, job_id)

    # # calls fortran with input file and generates .OUT files
    subprocess.call(
        [path + "boxm_orig", path, job_id],
    )

    cell_chem_ds = update_chem_ds(cell_chem_ds, job_id)

    return cell_chem_ds

def gen_boxm_orig_input(cell_chem_ds, job_id):
    """Generate the input file for the original box model."""

    # delete any existing input files
    if pathlib.Path(f"inputs/{job_id}/boxm_input.txt").exists():
            pathlib.Path(f"inputs/{job_id}/boxm_input.txt").unlink()

    # open file
    boxm_input = open(f"inputs/{job_id}/boxm_input.txt", "w")

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
    P = cell_chem_ds["air_pressure"].item()
    H2O = cell_chem_ds["H2O"].values[0]
    temp = cell_chem_ds["air_temperature"].values[0]

    boxm_input.write(f"{day}\n{month}\n{year}\n{level}\n{longbox}\n{latbox}\n{runtime}\n{M}\n{plevel}\n{H2O}\n{temp}\n")
    for s in ["NO2", "NO", "O3", "CO", "CH4", "HCHO", "CH3CHO", "CH3COCH3",
                        "C2H6", "C2H4", "C3H8", "C3H6", "C2H2", "NC4H10", "TBUT2ENE",
                        "BENZENE", "TOLUENE", "OXYL", "C5H8", "H2O2", "HNO3", "C2H5CHO",
                        "CH3OH", "MEK", "CH3OOH", "PAN", "MPAN"]:
        
        boxm_input.write(f"{cell_chem_ds["bg_chem"].sel(species=s).item()}\n")
        
    boxm_input.close()

def gen_zen_file(cell_chem_ds, job_id):
    """Generate the ZEN file for the original box model."""

    # delete any existing input files
    zen_file_path = pathlib.Path(f"inputs/{job_id}/zen.csv")
    if zen_file_path.exists():
        zen_file_path.unlink()

    # Extract the sza data and convert it to a DataFrame
    sza_data = cell_chem_ds["sza"].values
    sza_df = pd.DataFrame(sza_data, columns=["sza"])

    # Write the DataFrame to a CSV file
    sza_df.to_csv(zen_file_path, index=False, header=False)

def gen_emi_file(cell_chem_ds, job_id):
    """Generate the EMI file for the original box model."""

    # delete any existing input files
    emi_file_path = pathlib.Path(f"inputs/{job_id}/emi.csv")
    if emi_file_path.exists():
        emi_file_path.unlink()

    # Extract the emi data and convert it to a DataFrame
    emi_data = cell_chem_ds["emi"].values
    emi_df = pd.DataFrame(emi_data, columns=cell_chem_ds["emi_species"].values)

    # Write the DataFrame to a CSV file
    emi_df.to_csv(emi_file_path, index=False, header=False)

def latitude_to_latbox(latitude):
        # Map the latitude to the range 0-1
        normalized_latitude = (latitude + 87.5) / 180

        # Map the normalized latitude to the range 1-72
        latbox = normalized_latitude * 36 + 1

        # Round to the nearest integer and return
        return round(latbox)

def longitude_to_longbox(longitude):
        # Map the longitude to the range 0-1
        normalized_longitude = (longitude + 177.5) / 360

        # Map the normalized longitude to the range 1-144
        longbox = normalized_longitude * 72 + 1

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

def update_chem_ds(cell_chem_ds, job_id):
    sza_df = pd.read_csv(f"outputs/{job_id}/ZEN.OUT", header=0,
                        names=['TIME', 'ZEN'], dtype=np.float64)
        
    J_df = pd.read_csv(f"outputs/{job_id}/J.OUT", header=0,
                        names=['TIME', 'J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'J11', 'J12', 'J13','J14', 'J15', 'J16', 'J17', 'J18', 'J19', 'J20', 'J21', 'J22', 'J23', 'J24', 'J25', 'J26', 'J27', 'J28', 'J29', 'J30', 'J31', 'J32', 'J33', 'J34', 'J35', 'J36', 'J37', 'J38', 'J39', 'J40', 'J41', 'J42', 'J43', 'J44', 'J45', 'J46', 'J47', 'J48', 'J49', 'J50'], dtype=np.float64)

    DJ_df = pd.read_csv(f"outputs/{job_id}/DJ.OUT", header=0,
                            names=['TIME', 'DJ1', 'DJ2', 'DJ3', 'DJ4', 'DJ5', 'DJ6', 'DJ7', 'DJ8', 'DJ9', 'DJ10', 'DJ11', 'DJ12', 'DJ13','DJ14', 'DJ15', 'DJ16', 'DJ17', 'DJ18', 'DJ19', 'DJ20', 'DJ21', 'DJ22', 'DJ23', 'DJ24', 'DJ25', 'DJ26', 'DJ27', 'DJ28', 'DJ29', 'DJ30', 'DJ31', 'DJ32', 'DJ33', 'DJ34', 'DJ35', 'DJ36', 'DJ37', 'DJ38', 'DJ39', 'DJ40', 'DJ41', 'DJ42', 'DJ43', 'DJ44', 'DJ45', 'DJ46', 'DJ47', 'DJ48', 'DJ49', 'DJ50'], dtype=np.float64)

    RC_df = pd.read_csv(f"outputs/{job_id}/RC.OUT", header=0,
                        names=['TIME', 'RC1', 'RC2', 'RC3', 'RC4', 'RC5', 'RC6', 'RC7', 'RC8', 'RC9', 'RC10', 'RC11', 'RC12', 'RC13','RC14', 'RC15', 'RC16', 'RC17', 'RC18', 'RC19', 'RC20', 'RC21', 'RC22', 'RC23', 'RC24', 'RC25', 'RC26', 'RC27', 'RC28', 'RC29', 'RC30', 'RC31', 'RC32', 'RC33', 'RC34', 'RC35', 'RC36', 'RC37', 'RC38', 'RC39', 'RC40', 'RC41', 'RC42', 'RC43', 'RC44', 'RC45', 'RC46', 'RC47', 'RC48', 'RC49', 'RC50'], dtype=np.float64)

    # get species names
    header_names = ['TIME'] + list(cell_chem_ds["species"].values)

    Y_df = pd.read_csv(f"outputs/{job_id}/Y.OUT", header=0,
                            names=header_names, dtype=np.float64) 
    
    # # Update the chem_ds_stacked with the new data
    # Update zen data
    cell_chem_ds["sza_orig"] = (["time"], da.zeros((cell_chem_ds.sizes["time"])))
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
        
    cell_chem_ds["Y_orig"] = (["time", "species_out"], da.zeros((cell_chem_ds.sizes["time"], cell_chem_ds.sizes["species_out"])))
    for s, species_out in enumerate(cell_chem_ds["species_out"].values):
        cell_chem_ds["Y_orig"].loc[:, species_out] = Y_df[species_out].values
     
    return cell_chem_ds
    

def parse_args():
    parser = argparse.ArgumentParser(description="Overwrite parameters from command line")
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
    parser.add_argument("--lagrangian_tendency_of_air_pressure", type=float, help="Lagrangian tendency of air pressure in m/s")
    parser.add_argument("--species_in", type=str, help="Input species (comma-separated)")
    parser.add_argument("--species_out", type=str, help="Output species (comma-separated)")
    parser.add_argument("--gpat_path", type=str, help="GPAT directory")
    parser.add_argument("--job_id", type=str, help="Job ID")
    parser.add_argument("--run_gpat", action='store_true', help="Run the GPAT model")
       
    return parser.parse_args()

def update_sim_params_from_args(params, args):
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
    return cls(**dict_obj)