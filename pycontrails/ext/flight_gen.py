"""Sample trajectory points for formation flights."""

import numpy as np
import pandas as pd
import xarray as xr
import dask.array as da
import dataclasses
from pyproj import Geod
from pycontrails.core import models, Flight, Fleet, MetDataset, GeoVectorDataset
from pycontrails.models.ps_model import PSFlight 
from pycontrails.models.emissions import Emissions
from pycontrails.models.dry_advection import DryAdvection


# @dataclasses.dataclass
# class FlightGenParams:
#     """Parameters for formation flight trajectory plotting."""

#     t0_fl: pd.Timestamp = pd.Timestamp(2020, 1, 1, 12, 0)
#     rt_fl: pd.Timedelta = pd.Timedelta(minutes=30)
#     ts_fl: pd.Timedelta = pd.Timedelta(minutes=1)
#     ac_type: str = "A320"
#     fl0_speed: float = 100.0
#     fl0_heading: float = 0.0
#     fl0_coords0: tuple[float, float, float] = (0.0, 0.0, 0.0)
#     sep_dist: tuple[float, float, float] = (0.0, 0.0, 0.0)
#     n_ac: float = 2
    
class FlightGen:
    """Generate trajectories and estimate fuel burn and emissions from formation flights."""
    
    name = "traj_gen"
    long_name = "Formation Flight Trajectory Generator"
    

    def __init__(self, 
                 met: MetDataset,
                 fl_params: dict,
                 plume_params: dict,
                 chem_params: dict
                 ):
        self.met = met
        self.fl_params = fl_params
        self.plume_params = plume_params
        self.chem_params = chem_params

        self.flights = []

    def traj_gen(self) -> list[Flight]:
        """Generate formation flight trajectory points."""
        fl_params = self.fl_params
        flights = self.flights

        lon0, lat0, alt0 = fl_params["fl0_coords0"]
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
        fl0.attrs = {"flight_id": int(0), "aircraft_type": fl_params["ac_type"]}

        flights.append(fl0)

        fl = fl0

        # create follower flight trajectories
        for i in range(1, fl_params["n_ac"]):
            fl = fl.copy()
            

            # calculate new coords for follower flight
            lon_dx, lat_dx, _ = geod.fwd(lon0, lat0, heading, fl_params["sep_dist"][0])
            lon_dx_dy, lat_dx_dy, _ = geod.fwd(lon_dx, lat_dx, heading + 90, fl_params["sep_dist"][1])
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
            flights.append(fl)

            # Update starting coordinates for next flight
            lon0, lat0, alt0 = lon_dx_dy, lat_dx_dy, alt_dx_dy

        return flights
        
    def calc_fb_emissions(self) -> list[Flight]:
        """Calculate fuel burn and emissions for formation flights."""
        flights = self.flights
        met = self.met
        chem_params = self.chem_params

        ps_model = PSFlight()
        emi_model = Emissions()

        for i, fl in enumerate(flights):

            # downselect met data to the flight trajectory
            flights[i]["air_temperature"] = models.interpolate_met(met, fl, "air_temperature")
            flights[i]["specific_humidity"] = models.interpolate_met(met, fl, "specific_humidity")
            flights[i]["true_airspeed"] = fl.segment_groundspeed()
            
            # get ac performance data using Poll-Schumann Model
            flights[i] = ps_model.eval(flights[i])

            # get emissions data
            flights[i] = emi_model.eval(flights[i])

            # emission indices
            eis = {
            "co2_m": 3.16,
            "h2o_m": 1.23,
            "so2_m": 0.00084,
            "nvpm_m": flights[i]["nvpm_ei_m"],
            "nox_m": flights[i]["nox_ei"],
            "co_m": flights[i]["co_ei"],
            # hydrocarbon speciation
            "hcho_m": 0.12 * flights[i]["hc_ei"], # formaldehyde
            "ch3cho_m": 0.04 * flights[i]["hc_ei"], # acetaldehyde
            "c2h4_m": 0.15 * flights[i]["hc_ei"], # ethylene
            "c3h6_m": 0.04 * flights[i]["hc_ei"], # propene
            "c2h2_m": 0.04 * flights[i]["hc_ei"], # acetylene
            "benzene_m": 0.02 * flights[i]["hc_ei"] # benzene
            }

            # calculate emission mass per metre squared for each species
            for species, ei in eis.items():
                flights[i][species] = (ei * flights[i]["fuel_flow"] / flights[i]["true_airspeed"]) / chem_params["vres_chem"] / 1E+03    

        return flights
    
    def sim_plumes(self) -> list[pd.DataFrame]:
        """Simulate dry advection and dispersion for each flight in the formation."""
        met = self.met
        plume_params = self.plume_params
        flights = self.flights

        dry_adv = DryAdvection(met, plume_params)
        plumes = []

        fl_df = []
        pl_df = []
        
        for i, fl in enumerate(flights):
            
            pl = dry_adv.eval(fl)
            plumes.append(pl)

            # convert both flights and plumes to dataframes
            flights[i] = flights[i].dataframe
            plumes[i] = plumes[i].dataframe

            # calc plume heading
            plumes[i] = calc_heading(plumes[i])

            plumes[i]["flight_id"] = flights[i]["flight_id"][0]
            flights[i]["waypoint"] = flights[i].index

        # concatenate all flights and plumes into single dfs
        fl_df = pd.concat(flights)
        pl_df = pd.concat(plumes)

        # merge the two dataframes
        pl_df = pd.merge(
                    fl_df[['flight_id', 'waypoint', 'fuel_flow', 'true_airspeed', 'co2_m', 'h2o_m',
                    'so2_m', 'nox_m', 'co_m', 'hcho_m', 'ch3cho_m', 'c2h4_m', 'c3h6_m', 'c2h2_m', 'benzene_m', 'nvpm_m']],
                    pl_df[['flight_id', 'waypoint', 'time', 'longitude', 'latitude', 'level', 'width', 'heading']], on=['flight_id', 'waypoint']).sort_values(by=['time', 'flight_id', 'waypoint'])
    
        pl_df['sin_a'] = np.sin(np.radians(pl_df['heading']))
        pl_df['cos_a'] = np.cos(np.radians(pl_df['heading']))
        
        return pl_df
    

def calc_heading(pl_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate heading for each plume."""
    # Sort the dataframe by time and waypoint
    pl_df = pl_df.sort_values(by=["time", "waypoint"])

    # Group the dataframe by the timestep and apply the function
    pl_df['heading'] = pl_df.groupby('time').apply(calculate_heading_g).reset_index(drop=True)

    return pl_df

def calculate_heading_g(group):
    g = Geod(ellps="WGS84")

    startlat = group['latitude'].values[:-1]
    startlon = group['longitude'].values[:-1]
    endlat = group['latitude'].values[1:]
    endlon = group['longitude'].values[1:]
    az12, az21, dist = g.inv(startlon, startlat, endlon, endlat)

    heading = (90 - az12) % 360
    
    return pd.Series(np.concatenate([[np.nan], heading]), index=group.index)
        
       