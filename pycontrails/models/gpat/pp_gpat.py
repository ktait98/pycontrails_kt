import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import os
import re
from pycontrails.core import GeoVectorDataset
from pycontrails.models.gpat.gpat import mc_test, boxm_test
from dataclasses import dataclass, asdict, fields, is_dataclass
import cartopy.crs as ccrs
import cartopy.feature as cfeature

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

    def load_chem_ds(self, job_id, i_lat, i_lon, i_alt):
        outputs_dir = self.outputs_dir

        chem_ds = xr.open_dataset(outputs_dir + job_id + "/chem_" + job_id + ".nc")
        chem_ds = chem_ds.expand_dims(job_id=[job_id])
        chem_ds = chem_ds.assign_coords(species_out=chem_ds.attrs["species_out"])        
        print(f"Loaded chem ds for {job_id}")

        chem_ds = chem_ds.isel(job_id=0, latitude=i_lat, longitude=i_lon, level=i_alt)
        
        return chem_ds

    def load_chem_da(self, job_id, property):
        outputs_dir = self.outputs_dir

        chem_ds = xr.open_dataset(outputs_dir + job_id + "/chem_" + job_id + ".nc")
        chem_ds = chem_da.expand_dims(job_id=[job_id])
        
        chem_da = chem_ds[property]
        
        print(f"Loaded chem da for {job_id}, {property}")
        
        return chem_da


    ### Functions for calculating metrics
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

        kCH3O2_NO = 3.00E-12 * np.exp(280/T)**0.999
        kCH3O2_OH = 1.3E-10
        kCH3O2_HO2 = 4.10E-13 * np.exp(790/T)

        # should be αCH3O2 = (kCH3O2+NO × [NO] + kCH3O2+OH × [OH]) / 
        # ( kCH3O2+NO × [NO] + kCH3O2+OH × [OH] + kCH3O2+HO2 × [HO2])
        # but no reaction for CH3O2+OH in the mechanism
        alpha_CH3O2 = (kCH3O2_NO * NO + kCH3O2_OH * OH) / \
        (kCH3O2_NO * NO + kCH3O2_OH * OH + kCH3O2_HO2 * HO2)

        return alpha_CH3O2

    ### Functions for plotting
    # Function to create a bar chart
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
    def plot_scatter_plot(ax, data_list, x, y, labels, title, xlabel, ylabel, xticks=None, yticks=None, grid=True):
        for data, label in zip(data_list, labels):
            if isinstance(data, xr.DataArray):
                ax.scatter(data[x].values, data[y].values, label=label, color=lighter_green)
            else:
                ax.scatter(data[x], data[y], label=label, color=lighter_green)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)
        if grid:
            ax.grid(True)
        ax.legend()

    # Function to create a line plot with multiple datasets
    def plot_line_plot(ax, x, y, title, xlabel, ylabel, xticks=None, yticks=None, grid=True):
        ax.plot(x, y)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)
        if grid:
            ax.grid(True)
        ax.legend()

    # Function to create a spatial heatmap
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
