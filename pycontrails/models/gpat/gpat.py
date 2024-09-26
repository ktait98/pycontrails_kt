"""Gridded Plume Analysis Tool (GPAT). 

Simulate aircraft trajectories, estimate aircraft performance, fuel burn and emissions. 

Plot associated aircraft exhaust plumes, subject to Gaussian dispersion and advection. Aggregate plumes to an Eulerian grid for photochemical and microphysical processing."""

import os
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Geod
from pycontrails.core import Flight, GeoVectorDataset, MetDataArray, MetDataset, models
from pycontrails.core.models import Model, ModelParams
from dataclasses import dataclass
from typing import Tuple



### GPAT Model Parameters ###
@dataclass
class FlParams(ModelParams):
    """Default flight/fleet parameters."""
    t0_fl: pd.Timestamp  # flight start time
    rt_fl: pd.Timedelta  # flight run time
    ts_fl: pd.Timedelta  # flight time step
    ac_type: str  # aircraft type
    fl0_speed: float  # m/s
    fl0_heading: float  # deg
    fl0_coords0: Tuple[float, float, float]  # lat, lon, alt [deg, deg, m]
    sep_dist: Tuple[float, float, float]  # dx, dy, dz [m]
    n_ac: int  # number of aircraft

class PlumeParams(ModelParams):
    """Default plume dispersion parameters."""
    dt_integration: pd.Timedelta  # integration time step
    max_age: pd.Timedelta  # maximum age of the plume
    depth: float  # initial plume depth, [m]
    width: float  # initial plume width, [m]
    shear: float  # wind shear [1/s]

class MetParams(ModelParams):
    """Default meteorological parameters."""
    t0_met: pd.Timestamp  # meteorology start time
    rt_met: pd.Timedelta  # meteorology runtime
    ts_met: pd.Timedelta  # meteorology time step
    lat_bounds: Tuple[float, float]  # lat bounds [deg]
    lon_bounds: Tuple[float, float]  # lon bounds [deg]
    alt_bounds: Tuple[float, float]  # alt bounds [m]
    hres_met: float  # horizontal resolution [deg]
    vres_met: float  # vertical resolution [m]

class ChemParams(ModelParams):
    """Default chemistry parameters"""
    t0_chem: pd.Timestamp  # chemistry start time
    rt_chem: pd.Timedelta  # chemistry runtime
    ts_chem: pd.Timedelta  # chemistry time step
    lat_bounds: Tuple[float, float]  # lat bounds [deg]
    lon_bounds: Tuple[float, float]  # lon bounds [deg]
    alt_bounds: Tuple[float, float]  # alt bounds [m]
    hres_pl: float  # horizontal resolution of the plume, [deg]
    hres_chem: float  # horizontal resolution [deg]
    vres_chem: float  # vertical resolution [m]

class GPAT(Model):
    """Gridded Plume Analysis Tool (GPAT)."""
    def __init__(
            self, 
            fl_params: FlParams, 
            plume_params: PlumeParams, 
            met_params: MetParams,
            chem_params: ChemParams
            ):
        
        self.fl_params = fl_params
        self.plume_params = plume_params
        self.chem_params = chem_params
        self.met_params = met_params

    def eval(self):
        """Run the GPAT model."""

        # Generate formation flight trajectory points
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
        self.pl = self.sim_plumes()

        # Aggregate plumes to an Eulerian grid for photochemical and microphysical processing
        self.emi = self.plume_to_grid()

        # Run Contrail Model
        self.contrail = self.run_cc()

        # Run BOXM
        self.chem = self.run_boxm()


    def traj_gen(self) -> list[Flight]:
        """Generate formation flight trajectory points."""
        fl_params = self.fl_params
        flights = self.flights

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
                    (fl0["latitude"] > self.chem_params["lat_bounds"][0] + 0.05) & (fl0["latitude"] < self.chem_params["lat_bounds"][1] - 0.05) &
                    (fl0["longitude"] > self.chem_params["lon_bounds"][0] + 0.05) & (fl0["longitude"] < self.chem_params["lon_bounds"][1] - 0.05) &
                    (fl0["altitude"] > self.chem_params["alt_bounds"][0]) & (fl0["altitude"] < self.chem_params["alt_bounds"][1])

                )
        
        fl0 = fl0.filter(mask)
        flights.append(fl0)

        fl = fl0

        if fl_params["n_ac"] > 1:
            # create follower flight trajectories
            for i in range(1, fl_params["n_ac"]):
                fl = fl.copy()

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
                fl["latitude"] += dlat
                fl["longitude"] += dlon
                fl["altitude"] += dalt
                fl.attrs = {"flight_id": int(i), "aircraft_type": fl_params["ac_type"]}

                mask = (
                    (fl["latitude"] > self.chem_params["lat_bounds"][0] + 0.05) & (fl["latitude"] < self.chem_params["lat_bounds"][1] - 0.05) &
                    (fl["longitude"] > self.chem_params["lon_bounds"][0] + 0.05) & (fl["longitude"] < self.chem_params["lon_bounds"][1] - 0.05) &
                    (fl["altitude"] > self.chem_params["alt_bounds"][0]) & (fl["altitude"] < self.chem_params["alt_bounds"][1])
                )
                fl = fl.filter(mask)
                flights.append(fl)

                # Update starting coordinates for next flight
                lon0, lat0, alt0 = lon_dx_dy, lat_dx_dy, alt_dx_dy

        return flights
    
    def gen_met(self) -> MetDataset:
        """Generate meteorological data."""
        met_params = self.met_params
    
    def gen_bg_chem(self) -> xr.Dataset:
        """Generate background chemistry data."""
        chem_params = self.chem_params

    def ac_perf(self) -> list[Flight]:
        """Calculate aircraft performance using PS Model."""
        met = self.met
        flights = self.flights

    def emissions(self) -> list[Flight]:
        """Estimate emissions using Pycontrails Emissions Model."""
        met = self.met
        flights = self.flights

    def sim_plumes(self) -> list[pd.DataFrame]:
        """Simulate plume dispersion/advection using Pycontrails Dry Advection Model."""
        met = self.met
        flights = self.flights



    