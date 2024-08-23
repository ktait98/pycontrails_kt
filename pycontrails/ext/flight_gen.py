"""Sample trajectory points for formation flights."""

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.animation import FuncAnimation, PillowWriter
from pyproj import Geod
import os

from pycontrails.core import Flight, GeoVectorDataset, MetDataArray, MetDataset, models
from pycontrails.models.cocip import contrails_to_hi_res_grid
from pycontrails.models.dry_advection import DryAdvection
from pycontrails.models.emissions import Emissions
from pycontrails.models.ps_model import PSFlight
from pycontrails.physics import units


class FlightGen:
    """Generate trajectories and estimate fuel burn and emissions from formation flights."""

    name = "traj_gen"
    long_name = "Formation Flight Trajectory Generator"

    def __init__(self, met: MetDataset, fl_params: dict, plume_params: dict, chem_params: dict):
        self.met = met
        self.fl_params = fl_params
        self.plume_params = plume_params
        self.chem_params = chem_params

        self.flights = []

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
                flights.append(fl)

                # Update starting coordinates for next flight
                lon0, lat0, alt0 = lon_dx_dy, lat_dx_dy, alt_dx_dy

        return flights

    def calc_fb_emissions(self) -> list[Flight]:
        """Calculate fuel burn and emissions for formation flights."""
        flights = self.flights
        met = self.met
        plume_params = self.plume_params
        chem_params = self.chem_params

        ps_model = PSFlight()
        emi_model = Emissions()

        for i, fl in enumerate(flights):

            # downselect met data to the flight trajectory
            flights[i].downselect_met(met)
            # flights[i].intersect_met(met["specific_humidity"])
            flights[i]["air_temperature"] = models.interpolate_met(met, fl, "air_temperature")
            flights[i]["specific_humidity"] = models.interpolate_met(met, fl, "specific_humidity")
            flights[i]["true_airspeed"] = fl.segment_groundspeed()
                        
            # get ac performance data using Poll-Schumann Model
            flights[i] = ps_model.eval(flights[i])
            
            # get emissions data
            flights[i] = emi_model.eval(flights[i])

             # Iterate over the columns in the DataFrame
            for column in flights[i].dataframe.columns:
                # Replace NaN values in the column with the value from the previous row
                flights[i].dataframe[column] = flights[i].dataframe[column].fillna(method='ffill')

            # emission indices
            eis = {
                # primary combustion products
                "CO2": 3.16,
                "H2O": 1.23,
                "SO2": 0.00084,

                # secondary combustion products
                "nvPM": flights[i]["nvpm_ei_m"],
                "NO": 0.95 * flights[i]["nox_ei"],
                "NO2": 0.05 * flights[i]["nox_ei"],
                "CO": flights[i]["co_ei"],

                # hydrocarbon speciation
                "HCHO": 0.12 * flights[i]["hc_ei"],  # formaldehyde
                "CH3CHO": 0.04 * flights[i]["hc_ei"],  # acetaldehyde
                "C2H4": 0.15 * flights[i]["hc_ei"],  # ethylene
                "C3H6": 0.04 * flights[i]["hc_ei"],  # propene
                "C2H2": 0.04 * flights[i]["hc_ei"],  # acetylene
                "BENZENE": 0.02 * flights[i]["hc_ei"],  # benzene
            }

            # calculate emission mass per metre squared for each species
            for species, ei in eis.items():
               
                flights[i][species] = (
                    (ei * flights[i]["fuel_flow"] / flights[i]["true_airspeed"])
                    / chem_params["vres_chem"] / plume_params["width"]
                )

        return flights

    def sim_plumes(self) -> list[pd.DataFrame]:
        """Simulate dry advection and dispersion for each flight in the formation."""
        met = self.met
        plume_params = self.plume_params
        flights = self.flights

        dry_adv = DryAdvection(met, plume_params)

        plumes = []
        fl_df = []

        for i, fl in enumerate(flights):

            pl = dry_adv.eval(fl)
            plumes.append(pl)

            # convert both flights and plumes to dataframes
            flights[i] = flights[i].dataframe
            for column in flights[i].columns:
                # Replace NaN values in the column with the value from the previous row
                flights[i][column] = flights[i][column].fillna(method='ffill')
                
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

        pl_df["sin_a"] = np.sin(np.radians(pl_df["heading"]))
        pl_df["cos_a"] = np.cos(np.radians(pl_df["heading"]))

        pl_df["altitude"] = units.pl_to_m(pl_df["level"])

        self.plumes = pl_df

        return fl_df, pl_df

    def plume_to_grid(self, lats_pl, lons_pl, alts, times) -> MetDataset:
        """Convert plume data to a grid."""
        chem_params = self.chem_params
        pl_df = self.plumes

        # loop over time and plume property
        emi = xr.DataArray(
            np.zeros((len(lons_pl), len(lats_pl), len(alts), len(times), 9)),
        dims=["longitude", "latitude", "level", "time", "emi_species"],
        coords={
                "longitude": lons_pl,
                "latitude": lats_pl,
                "level": units.m_to_pl(alts),
                "time":  times,
                "emi_species": ["NO", "NO2", "CO", "HCHO", "CH3CHO", "C2H4", "C3H6", "C2H2", "BENZENE"],
            }
        )
        if self.fl_params["n_ac"] >= 1:
            
            for t, time in enumerate(pl_df["time"].unique()):
                print("Processing time: ", time)
                # create geovectordataset to store instantaneous plume data
                plume_time_data = GeoVectorDataset(data=pl_df.loc[pl_df["time"] == time])
                calc_continuous(plume_time_data)

                # define molar masses of species g/mol
                mm = [30.01, 46.01, 28.01, 30.03, 44.05, 28.05, 42.08, 26.04, 78.11]  # g/mol
                NA = 6.022e23  # Avogadro's number
                bbox = (
                    chem_params["lon_bounds"][0],
                    chem_params["lat_bounds"][0],
                    chem_params["lon_bounds"][1],
                    chem_params["lat_bounds"][1],
                    chem_params["alt_bounds"][0],
                    chem_params["alt_bounds"][1],
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
                        spatial_grid_res=chem_params["hres_pl"],
                    )

                    plume = (plume_property_data / 1E+03) * NA / mm[p] # kg/m^3 to molecules/cm^3 
                    # kg -> g (* 1E+03)
                    # m^3 -> cm^3 (/ 1E+06)

                    # Check if the 'plumes' directory exists, and create it if it does not
                    # if not os.path.exists("plumes"):
                    #     os.makedirs("plumes")
                    # np.savetxt("plumes/plume_data_" + repr(t) + "_" + repr(property) + ".csv", plume, delimiter=",") 

                    # find altitude index for flight level
                    emi.loc[:, :, units.m_to_pl(self.fl_params["fl0_coords0"][2]), time, property] = (
                        plume
                    )  # convert to molecules/cm^3

        return MetDataset(xr.Dataset({"emi": emi}))
    
    def anim_fl(self, fl_df: pd.DataFrame, pl_df: pd.DataFrame) -> None:
        """Animate formation flight trajectories and associated plume dispersion/advection."""

        fig1, ax1 = plt.subplots()

        scat_fl = ax1.scatter([], [], s=5, c="red", label="Flight path")
        scat_pl = ax1.scatter([], [], s=0.1, c="blue", label="Plume evolution")
        ax1.legend(loc="upper left")
        ax1.set_xlim([self.chem_params["lon_bounds"][0], self.chem_params["lon_bounds"][1]])
        ax1.set_ylim([self.chem_params["lat_bounds"][0], self.chem_params["lat_bounds"][1]])

        fl_frames = fl_df.groupby(fl_df["time"].dt.ceil(self.plume_params["dt_integration"]))
        pl_frames = pl_df.groupby(pl_df["time"].dt.ceil(self.plume_params["dt_integration"]))

        times = pl_frames.indices

        def animate(t):
            ax1.set_title(f"Time: {t}")

            try:
                group = fl_frames.get_group(t)
            except KeyError:
                offsets = [[None, None]]
            else:
                offsets = group[["longitude", "latitude"]]

            scat_fl.set_offsets(offsets)

            group = pl_frames.get_group(t)
            offsets = group[["longitude", "latitude"]]
            width = 10e-3 * group["width"]
            scat_pl.set_offsets(offsets)
            scat_pl.set_sizes(width)

            return scat_fl, scat_pl

        plt.close()
        anim = FuncAnimation(fig1, animate, frames=times)
        filename = pathlib.Path("plumes.gif")
        anim.save(filename, dpi=300, writer=PillowWriter(fps=8))


# functions used in flight gen
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
