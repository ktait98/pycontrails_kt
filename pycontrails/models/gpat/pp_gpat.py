import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
import os
import re
import subprocess
import pathlib
from pycontrails.core import GeoVectorDataset
from pycontrails.physics import units
from dataclasses import dataclass, asdict, fields, is_dataclass
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon
import matplotlib.patches as patches
from scipy.stats import multivariate_normal

class GPATPostProcessor:
    def __init__(self, outputs_dir, criteria):
        self.outputs_dir = outputs_dir
        self.criteria = criteria
        self.jobs_df = self.create_jobs_df()
        self.filtered_df = self.filter_jobs_df()
        self.job_ids = self.filtered_df.index.tolist()
        
    ### Functions for post-processing
    def create_jobs_df(self):
        outputs_dir = self.outputs_dir
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

    def filter_jobs_df(self):
        jobs_df = self.jobs_df
        criteria = self.criteria

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

    def load_fl_df(self, job_id):
        outputs_dir = self.outputs_dir

        fl_df = pd.read_pickle(outputs_dir + job_id + "/fl_" + job_id + ".pkl")
        fl_df['job_id'] = job_id  # Add job_id to the DataFrame

        print(f"Loaded flight data for {job_id}")

        fl_df = fl_df.set_index("job_id")

        return fl_df

    def load_pl_df(self, job_id):
        outputs_dir = self.outputs_dir

        pl_df = pd.read_pickle(outputs_dir + job_id + "/pl_" + job_id + ".pkl")
        pl_df['job_id'] = job_id  # Add job_id to the DataFrame

        pl_df = pl_df.set_index("job_id")  
        print(f"Loaded plume data for {job_id}")

        return pl_df

    def load_chem_ds(self, job_id, chunk_size=None):
        outputs_dir = self.outputs_dir

        # Define chunks if chunk_size is provided
        chunks = {'time': chunk_size} if chunk_size else None

        chem_ds = xr.open_dataset(outputs_dir + job_id + "/chem_" + job_id + ".nc", chunks=chunks)
        chem_ds = chem_ds.expand_dims(job_id=[job_id])
        chem_ds = chem_ds.assign_coords(species_out=chem_ds.attrs["species_out"])  
        chem_ds = chem_ds.isel(level=1)       
        print(f"Loaded chem ds for {job_id}")
        
        return chem_ds

    def load_chem_ds_cell(self, job_id, i_lat, i_lon):
        outputs_dir = self.outputs_dir

        chem_ds = xr.open_dataset(outputs_dir + job_id + "/chem_" + job_id + ".nc")
        chem_ds = chem_ds.expand_dims(job_id=[job_id])
        chem_ds = chem_ds.assign_coords(species_out=chem_ds.attrs["species_out"]) 
        chem_ds = chem_ds.isel(latitude=i_lat, longitude=i_lon)
        print(f"Loaded chem ds for {job_id}")
        
        return chem_ds

    def load_chem_ds_avg(self, job_id):
        outputs_dir = self.outputs_dir

        chem_ds = xr.open_dataset(outputs_dir + job_id + "/chem_" + job_id + ".nc")
        chem_ds = chem_ds.expand_dims(job_id=[job_id])
        chem_ds = chem_ds.assign_coords(species_out=chem_ds.attrs["species_out"])
        chem_ds = chem_ds.isel(altitude=1)
        chem_ds = chem_ds.mean(dim=["latitude", "longitude"])

    def load_chem_da(self, job_id, property):
        outputs_dir = self.outputs_dir

        chem_ds = xr.open_dataset(outputs_dir + job_id + "/chem_" + job_id + ".nc")
        chem_ds = chem_da.expand_dims(job_id=[job_id])
        
        chem_da = chem_ds[property]
        
        print(f"Loaded chem da for {job_id}, {property}")
        
        return chem_da


    ### Functions for calculating metrics
    # Function to calculate net ozone production rate
    def calc_NOPR(self, chem_ds):
        T = chem_ds["air_temperature"].isel(time=0)
        J1 = chem_ds["J"].isel(photol_params=0)
        H2O = chem_ds["H2O"]
        N2 = chem_ds["N2"]
        O2 = chem_ds["O2"]
        NO = chem_ds["Y"].sel(species_out="NO")
        OH = chem_ds["Y"].sel(species_out="OH")
        HO2 = chem_ds["Y"].sel(species_out="HO2")
        O3 = chem_ds["Y"].sel(species_out="O3")
        CH3O2 = chem_ds["Y"].sel(species_out="CH3O2")

        kO1D_H2O = 2.14E-10
        kO1D_N2 = 2.15E-11 * np.exp(110/T)
        kO1D_O2 = 3.2E-11 * np.exp(67/T)
        alpha_O1D = (kO1D_H2O * H2O) / (kO1D_N2 * N2 + kO1D_O2 * O2 + kO1D_H2O * H2O)

        kNO_HO2 = 8.1E-12
        kNO_CH3O2 = 2.3E-12 * np.exp(360/T)
        kO3_OH = 1.7E-12 * np.exp(940/T)
        kO3_HO2 =  2.03E-16 * ((T/400)**4.57) * np.exp(693/T)

        NOPR = NO * (kNO_HO2 * HO2 + kNO_CH3O2 * CH3O2) \
                - O3 * (kO3_HO2 * HO2 + kO3_OH * OH + alpha_O1D * J1)
        
        return NOPR

    # Function to calculate NOy
    def calc_NOy(self, chem_ds):
            NOy = chem_ds["Y"].sel(species_out="NO") + chem_ds["Y"].sel(species_out="NO2") + \
                chem_ds["Y"].sel(species_out="NO3") + chem_ds["Y"].sel(species_out="HNO3") + \
                chem_ds["Y"].sel(species_out="PAN")
            return NOy

    # Function to calculate NOz
    def calc_NOz(self, chem_ds):
        NOz = chem_ds["Y"].sel(species_out="HNO3") + chem_ds["Y"].sel(species_out="PAN") + \
            chem_ds["Y"].sel(species_out="NO3")

        return NOz

    # Function to calculate NO2t
    def calc_NO2t(self, chem_ds):
        NO = chem_ds["Y"].sel(species_out="NO")
        NO2 = chem_ds["Y"].sel(species_out="NO2")
        CO = chem_ds["Y"].sel(species_out="CO")
        CH4 = chem_ds["Y"].sel(species_out="CH4")
        OH = chem_ds["Y"].sel(species_out="OH")
        T = chem_ds["air_temperature"].isel(time=0)
        M = chem_ds["M"].isel(time=0)

        # Rate coefficients
        k1 = 5.4e-14 * (T / 298)**1.5 * np.exp(250 / T)
        k7 = 2.45e-12 * np.exp(-1775 / T)
        k0 = 2.5e-30 * (T / 300)**-4.4
        k_alpha = 1.6e-11 * (T / 300)**-1.7

        # Calculate k6[M]
        k6_M = (k0 * M) / (1 + (k0 * M) / k_alpha) * 0.6**(1 / (1 + np.log10(k0 * M / k_alpha)**2))

        NO2t = (k1 * CO + k7 * CH4) / k6_M

        return NO2t
    
    # Function to calculate O3_NOz
    def calc_O3_NOz(self, chem_ds, NOz):
        O3 = chem_ds["Y"].sel(species_out="O3")
        O3_NOz = O3 / NOz

        return O3_NOz
    
    # Function to calculate O3_NOy
    def calc_HCHO_NO2(self, chem_ds):
        HCHO = chem_ds["Y"].sel(species_out="HCHO")
        NO2 = chem_ds["Y"].sel(species_out="NO2")
        HCHO_NO2 = HCHO / NO2

        return HCHO_NO2
    
    # Function to calculate H2O2_HNO3
    def calc_H2O2_HNO3(self, chem_ds):
        H2O2 = chem_ds["Y"].sel(species_out="H2O2")
        HNO3 = chem_ds["Y"].sel(species_out="HNO3")
        H2O2_HNO3 = H2O2 / HNO3

        return H2O2_HNO3
    
    # Function to calculate CH3O2_alpha
    def calc_alpha_CH3O2(self, chem_ds):
        NO = chem_ds["Y"].sel(species_out="NO")
        OH = chem_ds["Y"].sel(species_out="OH")
        HO2 = chem_ds["Y"].sel(species_out="HO2")
        T = chem_ds["air_temperature"].isel(time=0)

        kCH3O2_NO = 2.30E-12 * np.exp(360/T)
        kCH3O2_OH = 3.7E-11 * np.exp(350/T) # 1.3E-10 # Assaf et al. 2016
        kCH3O2_HO2 = 3.8E-13 * np.exp(780/T)

        # should be αCH3O2 = (kCH3O2+NO × [NO] + kCH3O2+OH × [OH]) / 
        # ( kCH3O2+NO × [NO] + kCH3O2+OH × [OH] + kCH3O2+HO2 × [HO2])
        # but no reaction for CH3O2+OH in the mechanism
        alpha_CH3O2 = (kCH3O2_NO * NO + kCH3O2_OH * OH) / \
        (kCH3O2_NO * NO + kCH3O2_OH * OH + kCH3O2_HO2 * HO2)

        return alpha_CH3O2

    ### Functions to postprocess data
    # Function to calculate the mean and standard deviation of a variable
    def calc_cell(self, chem_ds, ilat, ilon):
        return chem_ds.isel(latitude=ilat, longitude=ilon)
    
    def calc_spatial_mean(self, da):
        return da.mean(dim=["latitude", "longitude"])
    
    def calc_temporal_mean(self, da):
        return da.mean(dim="time")

    def calc_daytime_mean(self, da, chem_ds):
        return da.where(chem_ds["sza"] < (np.pi/2)).mean(dim="time")
    
    def calc_nighttime_mean(self, da, chem_ds):
        return da.where(chem_ds["sza"] > (np.pi/2)).mean(dim="time")
    
    # Function to consider only plume affected cells
    def filter_plume_cells(self, chem_ds):
        # Check for non-zero emissions data
        plume_mask = chem_ds['emi'].sum(dim='emi_species') > 0

        return chem_ds.where(plume_mask)
    
    # Validation functions


    ### Functions for plotting
    # Function to create a bar chart
    @staticmethod
    def plot_bar_chart(ax, data, x, y, title, xlabel, ylabel, xticks=None, yticks=None, grid=True):
        ax.bar(data[x], data[y])
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)
        if grid:
            ax.grid(True)

    # Function to create a scatter plot with multiple datasets
    @staticmethod
    def plot_scatter_plot(ax, x, y, title, xlabel, ylabel, xticks=None, yticks=None, grid=True):
        ax.scatter(x, y)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)
        if grid:
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.legend()

    # Function to create a line plot with multiple datasets
    @staticmethod
    def plot_line_plot(ax, x, y, title, xlabel, ylabel, label=None, xticks=None, yticks=None, color='blue', grid=True, show_legend=True):
        ax.plot(x, y, color=color, label=label)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xlabel == "Time / days":
            # Format the x-ticks to show just the day number
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)
        if grid:
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        if label and show_legend:
            ax.legend(loc='upper left')

    # Function to create a spatial heatmap
    @staticmethod
    def plot_spatial_heatmap(ax, data, lon, lat, value, title, xlabel, ylabel, xticks=None, yticks=None, grid=True):
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAKES, alpha=0.5)
        ax.add_feature(cfeature.RIVERS)
        
        # Create a heatmap using pcolormesh
        heatmap = ax.pcolormesh(data[lon], data[lat], data[value], cmap='viridis', transform=ccrs.PlateCarree())
        plt.colorbar(heatmap, ax=ax, orientation='vertical', label=value)
        
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)
        if grid:
            ax.gridlines(draw_labels=True)

    # Data visualisation
    @staticmethod
    def plot_heatmap(params, fl_df, pl_df, chem_ds, **plot_params):
        fig1, ax1 = plt.subplots()
        ax1.set_xticks(np.arange(chem_ds["longitude"][0], chem_ds["longitude"][-1], 0.05))
        ax1.set_yticks(np.arange(chem_ds["latitude"][0], chem_ds["latitude"][-1], 0.05))

        print(f"ts: {plot_params['ts']}")
        print(f"time: {chem_ds["time"].values[plot_params["ts"]]}")

        # Plot the heatmap

        if plot_params['da'] == "emi":
            heatmap_data = (
                chem_ds[plot_params['da']].sel(emi_species=plot_params['property'], time=chem_ds["time"].values[plot_params["ts"]])
                .transpose("latitude", "longitude")
            )
        if plot_params['da'] == "Y":
            heatmap_data = (
                chem_ds[plot_params['da']].sel(species_out=plot_params['property'], time=chem_ds["time"].values[plot_params["ts"]])
                .transpose("latitude", "longitude")
            )

        # Shift the heatmap data up and right by half the resolution
        shifted_latitude = heatmap_data["latitude"] + params["hres_sim"] / 2
        shifted_longitude = heatmap_data["longitude"] + params["hres_sim"] / 2

        heatmap_data.coords["latitude"] = shifted_latitude
        heatmap_data.coords["longitude"] = shifted_longitude

        heatmap_data.plot(ax=ax1, cmap="summer")  # You can choose a colormap of your preference

        if plot_params['plot_fl'] == True:
            scat_fl = ax1.scatter(
                fl_df["longitude"].loc[fl_df["time"] == chem_ds["time"].values[plot_params["ts"]]],
                fl_df["latitude"].loc[fl_df["time"] == chem_ds["time"].values[plot_params["ts"]]],
                s=5,
                c="red",
                label="Flight path",
            )
        # Plot plume segments as polygons
        if plot_params['plot_pl'] == True:
            pl_df_time = pl_df[pl_df["time"] == chem_ds["time"].values[plot_params["ts"]]]
            for i, (index, row) in enumerate(pl_df_time.iterrows()):
                lon = row["longitude"]
                lat = row["latitude"]
                width = row["width"]
                heading = row["heading"]
                sigma_yy = row["sigma_yy"]

                # Skip if there is no next row
                if i >= len(pl_df_time) - 1:
                    continue

                next_row = pl_df_time.iloc[i + 1]
                next_lon = next_row["longitude"]
                next_lat = next_row["latitude"]
                next_width = next_row["width"]
                next_heading = next_row["heading"]


                dlon = units.m_to_longitude_distance(width * np.sin(np.radians(heading)) * 0.5, lat)
                next_dlon = units.m_to_longitude_distance(next_width * np.sin(np.radians(next_heading)) * 0.5, next_lat)
                dlat = units.m_to_latitude_distance(width * np.cos(np.radians(heading)) * 0.5)
                next_dlat = units.m_to_latitude_distance(next_width * np.cos(np.radians(next_heading)) * 0.5)

                corners = [
                    (lon - dlon, lat - dlat),
                    (next_lon + next_dlon, next_lat - next_dlat),
                    (next_lon - next_dlon, next_lat + next_dlat),
                    (lon + dlon, lat + dlat),

                ]

                polygon = Polygon(corners)
                patch = patches.Polygon(np.array(polygon.exterior.coords), closed=True, color='blue', facecolor='none', alpha=0.1)
                ax1.add_patch(patch)

                # # Generate Gaussian distribution
                # x = np.linspace(lon - 3 * sigma_yy, lon + 3 * sigma_yy, 100)
                # y = np.linspace(lat - 3 * sigma_yy, lat + 3 * sigma_yy, 100)
                # X, Y = np.meshgrid(x, y)
                # pos = np.dstack((X, Y))
                # rv = multivariate_normal([lon, lat], [[sigma_yy, 0], [0, sigma_yy]])
                # Z = rv.pdf(pos)

                # # Overlay Gaussian heatmap
                # ax1.contourf(X, Y, Z, levels=10, cmap="Blues", alpha=0.5)

        ax1.legend(loc="upper left")
        ax1.set_xlim([params["lon_bounds"][0], params["lon_bounds"][1]])
        ax1.set_ylim([params["lat_bounds"][0], params["lat_bounds"][1]])
        plt.grid()
        plt.show()

    @staticmethod
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


def mc_test(params, fl_df, pl_df, chem_ds):
    """Check if mass is conserved in the box model."""

    # Initialize the dictionary
    vecmass = {emi_species: [] for emi_species in chem_ds["emi_species"].values.tolist()}
    gridmass = {emi_species: [] for emi_species in chem_ds["emi_species"].values.tolist()}
    mc = {emi_species: [] for emi_species in chem_ds["emi_species"].values.tolist()}

    # Constants
    mm = [30.01, 46.01, 28.01, 30.03, 44.05, 28.05, 42.08, 26.04, 78.11]  # g/mol
    NA = 6.022e23  # Avogadro's number

    for s, emi_species in enumerate(chem_ds["emi_species"].values):
        
        max_fl_time = fl_df["time"].max()

        for ts, time in enumerate(pl_df["time"].unique()[:-1]):
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

            if time <= max_fl_time:
                # Accumulate vector mass for all flights
                vector_mass = fl_snapshot[emi_species]

                total_vector_mass += vector_mass.sum()

            # Grab plume mass from grid data
            grid_concs = chem_ds["emi"].sel(emi_species=emi_species, time=time)

            if (grid_concs == 0).all():
                pass
            else:
                # Compute the boolean indexer first
                grid_concs_over_zero = (grid_concs > 0)
                grid_concs_over_zero = grid_concs.where(grid_concs_over_zero, drop=True)
                               
                grid_mass = grid_concs_over_zero \
                    * chem_ds["M"].sel(time=time) \
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
    vecmass = pd.DataFrame(vecmass, index=pl_df["time"].unique()[:-1], columns=chem_ds["emi_species"].values.tolist())
    gridmass = pd.DataFrame(gridmass, index=pl_df["time"].unique()[:-1], columns=chem_ds["emi_species"].values.tolist())
    mc = pd.DataFrame(mc, index=pl_df["time"].unique()[:-1], columns=chem_ds["emi_species"].values.tolist())

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
    
def r_sq(y_true, y_pred):
    """Calculate the R-squared value for a model."""
    y_true_mean = np.mean(y_true)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true_mean) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    return r2