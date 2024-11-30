import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import os
import re
from pycontrails.core import GeoVectorDataset
from pycontrails.models.gpat.gpat import create_jobs_df, filter_jobs_df, load_fl_df, load_pl_df, load_chem_ds, load_chem_da, mc_test, boxm_test
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def main():
    outputs_dir = f"{os.getcwd()}/outputs/"

    # Filter criteria
    criteria = {
            # "n_ac": 3,
            # "rt_fl": (pd.Timedelta(minutes=30), pd.Timedelta(hours=2)),
            # "date_created": (pd.Timestamp("2024-11-16"), pd.Timestamp("2024-11-17")),
            "job_id": ['id_vs_plume_NA_1_1000_0_2_0.01_0.05']

        }

    # Load data
    jobs_df = create_jobs_df(outputs_dir)

    # Filter data
    filtered_df = filter_jobs_df(jobs_df, criteria)

    job_ids = filtered_df.index.values
    
    fl_df = load_fl_df(job_ids, outputs_dir)
    pl_df = load_pl_df(job_ids, outputs_dir)
    chem_ds = load_chem_ds(job_ids, outputs_dir)

    # BG run plots
    # Example NA background run?
    fig1, ax1 = plt.subplots(1, 1, figsize=(15, 10))
    plot_line_plot(ax1, fl_df, x='time', y='altitude', title='', xlabel='Time', ylabel='Altitude')

    # Show scatter plot of metrics across all months (one subplot for each hotspot)

    # Show scatter plot of metrics across all hotspots, averaged over all four months

    # Adjust layout
    plt.tight_layout()

    # Show the plots
    plt.show()

    # Save the plots
    def save_plot(filename, format='png'):
        plt.savefig(f'{filename}.{format}', format=format)

    # Example of saving a plot
    save_plot('/user/home/kt16229/work/pycontrails_kt/pycontrails/models/gpat/outputs/plots/multiple_plots', format='png')


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
        ax.scatter(data[x], data[y], label=label)
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
def plot_line_plot(ax, data_list, x, y, labels, title, xlabel, ylabel, xticks=None, yticks=None, grid=True):
    for data, label in zip(data_list, labels):
        ax.plot(data[x], data[y], label=label)
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