"""Gridded Plume Analysis Tool (GPAT). 

Simulate aircraft trajectories, estimate aircraft performance, fuel burn and emissions. 

Plot associated aircraft exhaust plumes, subject to Gaussian dispersion and advection. Aggregate plumes to an Eulerian grid for photochemical and microphysical processing."""

import os
import numpy as np
import pandas as pd
import xarray as xr
import dask.array as da
from pyproj import Geod
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import subprocess
from pycontrails.core import Flight, GeoVectorDataset, MetDataArray, MetDataset, models
from pycontrails.core.models import Model, ModelParams
from pycontrails.models.emissions import Emissions
from pycontrails.models.ps_model import PSFlight
from pycontrails.models.dry_advection import DryAdvection
from pycontrails.models.cocip import contrails_to_hi_res_grid
from pycontrails.physics import geo, thermo, units, constants
from dataclasses import dataclass
from typing import Tuple
import pathlib

### GPAT Model Parameters ###
@dataclass
class FlParams():
    """Default flight/fleet parameters."""
    t0_fl: pd.Timestamp = pd.to_datetime("2022-01-20 13:00:00") # flight start time
    rt_fl: pd.Timedelta = pd.Timedelta(minutes=60) # flight run time
    ts_fl: pd.Timedelta = pd.Timedelta(minutes=2) # flight time step
    ac_type: str = "A320" # aircraft type
    fl0_speed: float = 100.0 # m/s
    fl0_heading: float = 0.0 # deg
    fl0_coords0: Tuple[float, float, float] = (0.1, 0.125, 12500) # lat, lon, alt [deg, deg, m]
    sep_dist: Tuple[float, float, float] = (5000, 2000, 0) # dx, dy, dz [m]
    n_ac: int = 1 # number of aircraft

class PlumeParams(ModelParams):
    """Default plume dispersion parameters."""
    dt_integration: pd.Timedelta = pd.Timedelta(minutes=2) # integration time step
    max_age: pd.Timedelta = pd.Timedelta(hours=2) # maximum age of the plume
    depth: float = 50.0 # initial plume depth, [m]
    width: float = 50.0 # initial plume width, [m]
    shear: float = 0.01 # wind shear [1/s]
    hres_pl: float = 0.01 # horizontal resolution of the plume, [deg]
    vres_pl: float = 500 # vertical resolution of the plume [m]

class SimParams(ModelParams):
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
    species_out: np.array = np.array(["O3", "NO2", "NO", 
                                    "NO3", "N2O5", "HNO3", 
                                    "HONO", "HO2", "OH", 
                                    "H2O2", "H2O", "CO", 
                                    "CH4", "C2H6", "C3H8", 
                                    "C2H4", "C3H6"])

class GPAT(Model):
    """Gridded Plume Analysis Tool (GPAT).
    
    Simulate aircraft trajectories, estimate aircraft performance, fuel burn and emissions. Then aggregates emissions, bg chemistry and meteorology to an Eulerian grid for photochemical and microphysical processing.
    
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

        # Set the model parameters
        self.fl_params = fl_params
        self.plume_params = plume_params
        self.sim_params = sim_params

        # Generate the grid
        self.lats_pl = np.arange(
            sim_params["lat_bounds"][0], sim_params["lat_bounds"][1] + plume_params["hres_pl"], plume_params["hres_pl"]
        )
        self.lons_pl = np.arange(
            sim_params["lon_bounds"][0], sim_params["lon_bounds"][1] + plume_params["hres_pl"], plume_params["hres_pl"]
        )
        self.lats = np.arange(
            sim_params["lat_bounds"][0], sim_params["lat_bounds"][1] + sim_params["hres_sim"], sim_params["hres_sim"]
        )
        self.lons = np.arange(
            sim_params["lon_bounds"][0], sim_params["lon_bounds"][1] + sim_params["hres_sim"], sim_params["hres_sim"]
        )
        self.alts = np.arange(
            sim_params["alt_bounds"][0], sim_params["alt_bounds"][1] + sim_params["vres_sim"], sim_params["vres_sim"]
        )
        self.levels = units.m_to_pl(self.alts)

        self.times = pd.date_range(
            start=sim_params["t0_sim"],
            end=sim_params["t0_sim"] + sim_params["rt_sim"],
            freq=sim_params["ts_sim"],
        )

        # Set the path to the model and files
        self.path = "/home/ktait98/pycontrails_kt/pycontrails/models/gpat/" 
        #os.environ['PYCONTRAILSDIR'] + "/models/boxmodel/"
        self.inputs = "/home/ktait98/pycontrails_kt/pycontrails/models/gpat/inputs/" 
        #os.environ['PYCONTRAILSDIR'] + "/models/files/"
        self.outputs = "/home/ktait98/pycontrails_kt/pycontrails/models/gpat/outputs/" 
        #os.environ['PYCONTRAILSDIR'] + "/models/files/"

    def eval(self):
        """Run the GPAT model."""

        # Generate formation flight trajectory points
        self.fl = self.traj_gen()

        # Generate meteorological data
        self.met = self.gen_met()
        
        # # Generate background chemistry data
        self.bg_chem = self.gen_bg_chem()

        # Calculate aircraft performance using PS Model
        self.fl = self.ac_perf()

        # Estimate emissions using Pycontrails Emissions Model
        self.fl = self.emissions()

        # Simulate plume dispersion/advection using Pycontrails Dry Advection Model
        self.pl = self.sim_plumes()

        # Aggregate plumes to an Eulerian grid for photochemical and microphysical processing
        self.emi = self.plume_to_grid()

        # # Run Contrail Model
        # # self.contrail = self.run_cc()

        # Run BOXM
        self.chem = self.run_boxm()

    # Model methods
    def traj_gen(self) -> list[Flight]:
        """Generate formation flight trajectory points."""
        fl_params = self.fl_params
        fl = []

        lat0, lon0, alt0 = fl_params["fl0_coords0"]
        heading = fl_params["fl0_heading"]
        dist = fl_params["fl0_speed"] * fl_params["rt_fl"].total_seconds()

        # calculate the final coordinates
        geod = Geod(ellps="WGS84")
        lon1, lat1, _ = geod.fwd(lon0, lat0, heading, dist)

        # create flight object for leader flight and resample points according to ts_fl
        df = pd.DataFrame()
        df["longitude"] = [lon0, lon1]
        df["latitude"] = [lat0, lat1]
        df["altitude"] = [alt0, alt0]
        df["time"] = [fl_params["t0_fl"], (fl_params["t0_fl"] + fl_params["rt_fl"])]

        ts_fl_min = int(fl_params["ts_fl"].total_seconds() / 60)

        fl0 = Flight(df).resample_and_fill(freq=f"{ts_fl_min}min")
        fl0.attrs = {"flight_id": 0, "aircraft_type": fl_params["ac_type"]}
        mask = (
                    (fl0["latitude"] > self.sim_params["lat_bounds"][0] + 0.05) & (fl0["latitude"] < self.sim_params["lat_bounds"][1] - 0.05) &
                    (fl0["longitude"] > self.sim_params["lon_bounds"][0] + 0.05) & (fl0["longitude"] < self.sim_params["lon_bounds"][1] - 0.05) &
                    (fl0["altitude"] > self.sim_params["alt_bounds"][0]) & (fl0["altitude"] < self.sim_params["alt_bounds"][1])

                )
        
        fl0 = fl0.filter(mask)
        fl.append(fl0)

        fli = fl0

        if fl_params["n_ac"] > 1:
            # create follower flight trajectories
            for i in range(1, fl_params["n_ac"]):
                fli = fli.copy()

                # calculate new coords for follower flight
                lon_dx, lat_dx, _ = geod.fwd(lon0, lat0, heading, fl_params["sep_dist"][0])
                lon_dx_dy, lat_dx_dy, _ = geod.fwd(
                    lon_dx, lat_dx, heading + 90, fl_params["sep_dist"][1]
                )
                alt_dx_dy = alt0 + fl_params["sep_dist"][2]

                # Calculate the differences in lat, lon, alt
                dlat = lat_dx_dy - lat0
                dlon = lon_dx_dy - lon0
                dalt = alt_dx_dy - alt0

                # Update the latitude and longitude of each point in the flight path
                fli["latitude"] += dlat
                fli["longitude"] += dlon
                fli["altitude"] += dalt
                fli.attrs = {"flight_id": int(i), "aircraft_type": fl_params["ac_type"]}

                mask = (
                    (fli["latitude"] > self.sim_params["lat_bounds"][0] + 0.05) & (fl["latitude"] < self.sim_params["lat_bounds"][1] - 0.05) &
                    (fli["longitude"] > self.sim_params["lon_bounds"][0] + 0.05) & (fl["longitude"] < self.sim_params["lon_bounds"][1] - 0.05) &
                    (fli["altitude"] > self.sim_params["alt_bounds"][0]) & (fl["altitude"] < self.sim_params["alt_bounds"][1])
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
                "eastward_wind": (("latitude", "longitude", "level", "time"), da.full((len(self.lats), len(self.lons), len(self.alts), len(self.times)), sim_params["eastward_wind"])),
                "northward_wind": (("latitude", "longitude", "level", "time"), da.full((len(self.lats), len(self.lons), len(self.alts), len(self.times)), sim_params["northward_wind"])),
                "lagrangian_tendency_of_air_pressure": (("latitude", "longitude", "level", "time"), da.full((len(self.lats), len(self.lons), len(self.alts), len(self.times)), sim_params["lagrangian_tendency_of_air_pressure"])),
                "air_temperature": (("latitude", "longitude", "level", "time"), da.zeros((len(self.lats), len(self.lons), len(self.alts), len(self.times)))),
            },

            coords={
                "longitude": self.lons, "latitude": self.lats, "level": units.m_to_pl(self.alts), "time": self.times
            }
        )

        met = MetDataset(met)

        month = self.times[0].month

        air_temperature = xr.open_dataarray(
            self.inputs + "air_temperature.nc"
        ).sel(month=month - 1).interp(
            longitude=self.lons, latitude=self.lats, level=self.levels,
            method="linear").broadcast_like(met.data)
        #os.environ['PYCONTRAILSDIR'] + "/models/boxmodel/air_temperature.nc"
        
        h2o_concs = xr.open_dataarray(
            self.inputs + "h2o_concs.nc"
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
        sim_params = self.sim_params

        month = self.times[0].month

        bg_chem = xr.open_dataarray(
            self.inputs + "species.nc"
        ).sel(month=month - 1)

        for s in [1, 2, 3, 5, 7, 9, 10, 13, 15, 16, 17, 18, 19, 20, 22, 24, 26, 27, 29, 31, 33, 35, 36, 37, 38, 40, 41, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 62, 63, 65, 66, 68, 69, 70, 72, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 200, 201, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219]:
            bg_chem[:, :, :, s-1] = 0

        bg_chem = bg_chem * 1e09  # convert mixing ratio to ppb

        # downselect and interpolate bg_chem to the simulation grid
        bg_chem = bg_chem.interp(
            longitude=self.lons, latitude=self.lats, level=self.levels
    )

        return bg_chem

    def ac_perf(self) -> list[Flight]:
        """Calculate aircraft performance using PS Model."""
        met = self.met
        fl = self.fl

        ps_model = PSFlight()

        for i, fli in enumerate(fl):
            # downselect met data to the flight trajectory
            fl[i].downselect_met(met)
            fl[i]["air_temperature"] = models.interpolate_met(met, fli, "air_temperature")
            fl[i]["specific_humidity"] = models.interpolate_met(met, fli, "specific_humidity")
            fl[i]["true_airspeed"] = fli.segment_groundspeed()
        
            # get ac performance data using Poll-Schumann Model
            fl[i] = ps_model.eval(fl[i])

        return fl

    def emissions(self) -> list[Flight]:
        """Estimate emissions using Pycontrails Emissions Model."""
        plume_params = self.plume_params
        met = self.met
        fl = self.fl

        emi_model = Emissions()

        for i, fli in enumerate(fl):

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
                    ei * fl[i]["fuel_flow"] * plume_params["dt_integration"].seconds
                )

        return fl

    def sim_plumes(self) -> list[pd.DataFrame]:
        """Simulate plume dispersion/advection using Pycontrails Dry Advection Model."""
        plume_params = self.plume_params
        met = self.met
        fl = self.fl

        # Create a new dictionary excluding hres_pl and vres_pl to input to dry advection model
        filtered_plume_params = {key: value for key, value in plume_params.items() if key not in {"hres_pl", "vres_pl"}}

        dry_adv = DryAdvection(met, **filtered_plume_params)

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
        pl = pd.merge(
            fl_df[
                [
                    "flight_id",
                    "waypoint",
                    "fuel_flow",
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
                    "longitude",
                    "latitude",
                    "level",
                    "width",
                    "heading",
                ]
            ],
            on=["flight_id", "waypoint"],
        ).sort_values(by=["time", "flight_id", "waypoint"])

        pl["sin_a"] = np.sin(np.radians(pl["heading"]))
        pl["cos_a"] = np.cos(np.radians(pl["heading"]))
        pl["altitude"] = units.pl_to_m(pl["level"])

        return pl

    def plume_to_grid(self) -> MetDataset:
        """Aggregate plumes to an Eulerian grid for photochemical and microphysical processing."""
        plume_params = self.plume_params
        pl = self.pl

        # loop over time and plume property
        emi = xr.DataArray(
            np.zeros((len(self.lons_pl), len(self.lats_pl), len(self.alts), len(self.times), 9)),
        dims=["longitude", "latitude", "level", "time", "emi_species"],
        coords={
                "longitude": self.lons_pl,
                "latitude": self.lats_pl,
                "level": units.m_to_pl(self.alts),
                "time":  self.times,
                "emi_species": ["NO", "NO2", "CO", "HCHO", "CH3CHO", "C2H4", "C3H6", "C2H2", "BENZENE"],
            }
        )

        if self.fl_params["n_ac"] >= 1:
            
            for t, time in enumerate(pl["time"].unique()):
                print("Processing time: ", time)
                # create geovectordataset to store instantaneous plume data
                plume_time_data = GeoVectorDataset(data=pl.loc[pl["time"] == time])
                calc_continuous(plume_time_data)

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

                for p, property in enumerate(["NO", "NO2", "CO"]):
                    # ,
                    # 'hcho_m',
                    # 'ch3cho_m',
                    # 'c2h4_m',
                    # 'c3h6_m',
                    # 'c2h2_m',
                    # 'benzene_m']):

                    # call contrails_to_hi_res_grid
                    plume_property_data = contrails_to_hi_res_grid(
                        time=time,
                        contrails_t=plume_time_data,
                        var_name=property,
                        spatial_bbox=bbox,
                        spatial_grid_res=plume_params["hres_pl"],
                    )

                    # convert mass to density [kg/m^3]
                    density = plume_property_data / (plume_params["vres_pl"] \
                    * units.latitude_distance_to_m(plume_params["hres_pl"]) \
                    * units.longitude_distance_to_m(plume_params["hres_pl"], (self.lats[0] + self.lats[-1]) / 2))

                    plume = (density / 1E+03) * NA / mm[p] # [kg/m^3] to [molecules/cm^3] 
                    # kg -> g (* 1E+03)
                    # m^3 -> cm^3 (/ 1E+06)

                    # find altitude index for flight level
                    emi.loc[:, :, units.m_to_pl(self.fl_params["fl0_coords0"][2]), time, property] = (
                        plume
                    )  # convert to molecules/cm^3

        emi = MetDataset(xr.Dataset({"emi": emi}))

        return emi

    def run_cc(self) -> xr.Dataset:
        """Run Contrail Model."""
        met = self.met
        emi = self.emi

        return contrail
    
    def run_boxm(self) -> xr.Dataset:
        """Run BOXM."""
        # Initialize the box model dataset
        self.init_boxm_ds()
        # Convert the datasets to netCDF
        self.to_netcdf()
        # Run the box model
        self.do_boxm()
        # Unstack the box model dataset
        self.unstack()

        chem = self.boxm_ds_unstacked

        # Save the box model dataset to a netCDF file
        chem.to_netcdf(self.outputs + "chem.nc")

        return chem
    
    # Methods for running the box model
    def init_boxm_ds(self):
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
            dts=self.sim_params["ts_sim"].total_seconds(),
            species_out=self.sim_params["species_out"],
            )

        self.boxm_ds["J"] = (["time", "level", "longitude", "latitude", "photol_params"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], 5)))

        self.boxm_ds["DJ"] = (["time", "level", "longitude", "latitude", "photol_coeffs"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], 5)))

        self.boxm_ds["RC"] = (["time", "level", "longitude", "latitude", "therm_coeffs"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], 5)))

        self.boxm_ds["Y"] = (["time", "level", "longitude", "latitude", "species_out"], da.zeros((self.boxm_ds.dims["time"], self.boxm_ds.dims["level"], self.boxm_ds.dims["longitude"], self.boxm_ds.dims["latitude"], len(self.sim_params["species_out"]))))

    def to_netcdf(self):
        """Convert the met, bg_chem, and emi datasets to boxm_ds.nc for use in the box model."""

        # stack datasets to get cell index for fortran
        self.boxm_ds_stacked = self.boxm_ds.stack(
            {"cell": ["level", "longitude", "latitude"]}
        )
        print(self.boxm_ds_stacked.indexes["cell"])

        self.boxm_ds_stacked = self.boxm_ds_stacked.reset_index("cell")

        # Delete any existing netCDF files
        if pathlib.Path(self.inputs + "boxm_ds.nc").exists():
            print("deleting boxm_ds.nc")
            pathlib.Path(self.inputs + "boxm_ds.nc").unlink()

        # Convert DataFrames to Datasets and write to netCDF
        print("WARNING: " + self.inputs + "boxm_ds.nc")
        self.boxm_ds_stacked.to_netcdf(self.inputs + "boxm_ds.nc", mode="w")

        print(self.boxm_ds_stacked)

    def do_boxm(self):
        """Run the box model in fortran using subprocess."""
        
        subprocess.call(
            [self.path + "boxm"]
        )
        
        # open nc file
        self.boxm_ds = xr.open_dataset(self.inputs + "boxm_ds.nc")
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

        filename = pathlib.Path(self.outputs + var1 + "_" + var2 + ".gif")

        anim.save(filename, dpi=300, writer=PillowWriter(fps=8))


# Functions used in GPAT Model

def calc_heading(pl_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate heading for each plume."""
    # Sort the dataframe by time and waypoint
    pl_df = pl_df.sort_values(by=["time", "waypoint"])

    # Group the dataframe by the timestep and apply the function
    pl_df["heading"] = pl_df.groupby("time").apply(calculate_heading_g).reset_index(drop=True)

    return pl_df

def calculate_heading_g(group):
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
    """Calculate szas for each cell at all timesteps."""
    sza = np.zeros((len(latitudes), len(longitudes), len(timesteps)))

    for lon, lonval in enumerate(longitudes):
        for lat, latval in enumerate(latitudes):

            theta_rad = geo.orbital_position(timesteps)

            sza[lat, lon, :] = np.arccos(
                geo.cosine_solar_zenith_angle(lonval, latval, timesteps, theta_rad)
            )
    return sza

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