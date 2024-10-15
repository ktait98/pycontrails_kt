import numpy as np
import numpy.typing as npt
import pandas as pd
import xarray as xr
import warnings
import matplotlib.pyplot as plt
from scipy.stats import norm
from shapely.geometry import Point, Polygon, box 
from shapely import intersection
from pycontrails.core.met import MetDataArray, MetDataset
from pycontrails.core.vector import GeoVectorDataset, vector_to_lon_lat_grid
from pycontrails.models.cocip.contrail_properties import contrail_edges, plume_mass_per_distance

from pycontrails.physics import geo, thermo, units
from pycontrails.utils import dependencies

def plume_to_grid(
    time: pd.Timestamp | np.datetime64,
    plumes_t: GeoVectorDataset,
    *,
    var_name: str,
    spatial_bbox: tuple[float, float, float, float] = (-180.0, -90.0, 180.0, 90.0),
    spatial_grid_res: float = 0.05,
) -> xr.DataArray:
    r"""
    Aggregate plume segments to a high-resolution longitude-latitude grid.

    Parameters
    ----------
    time : pd.Timestamp | np.datetime64
        UTC time of interest.
    plumes_t : GeoVectorDataset
        All plume waypoint outputs at `time`.
    var_name : str
        Plume property for aggregation, where `var_name` must be included in `plume_segment`.
        For example, `tau_contrail`, `rf_sw`, `rf_lw`, and `rf_net`
    spatial_bbox : tuple[float, float, float, float]
        Spatial bounding box, `(lon_min, lat_min, lon_max, lat_max)`, [:math:`\deg`]
    spatial_grid_res : float
        Spatial grid resolution, [:math:`\deg`]

    Returns
    -------
    xr.DataArray
        Plume segments and their properties aggregated to a longitude-latitude grid.
    """
    # Ensure the required columns are included in `plumes_t`
    cols_req = [
        "flight_id",
        "waypoint",
        "longitude",
        "latitude",
        "altitude",
        "time",
        "sin_a",
        "cos_a",
        "width",
        "sigma_yy",
        var_name,
    ]
    plumes_t.ensure_vars(cols_req)

    # Ensure that the times in `plumes_t` are the same.
    is_in_time = plumes_t["time"] == time
    if not np.all(is_in_time):
        warnings.warn(
            f"Plume segments have inconsistent times. Waypoints that are not in {time} are removed."
        )
        plumes_t = plumes_t.filter(is_in_time)

    main_grid = _initialise_longitude_latitude_grid(spatial_bbox, spatial_grid_res)

    # Plume head and tails: continuous segments only
    heads_t = plumes_t.dataframe
    heads_t = heads_t.sort_values(["flight_id", "waypoint"])
    tails_t = heads_t.shift(periods=-1)

    is_continuous = heads_t["continuous"]
    heads_t = heads_t[is_continuous].copy()
    tails_t = tails_t[is_continuous].copy()
    tails_t["waypoint"] = tails_t["waypoint"].astype("int")

    heads_t = heads_t.set_index(["flight_id", "waypoint"], drop=False)
    tails_t.index = heads_t.index

    # Aggregate plume segments to a high resolution longitude-latitude grid
    try:
        from tqdm.auto import tqdm
    except ModuleNotFoundError as exc:
        dependencies.raise_module_not_found_error(
            name="plume_to_grid function",
            package_name="tqdm",
            module_not_found_error=exc,
        )

    for i in tqdm(heads_t.index):
        plume_segment = GeoVectorDataset(
            pd.concat([heads_t[cols_req].loc[i], tails_t[cols_req].loc[i]], axis=1).T, copy=True
        )

        segment_grid = segment_property_to_hi_res_grid(
            plume_segment, 
            var_name=var_name, 
            spatial_grid_res=spatial_grid_res,
            n_slices=5
        )
        main_grid = _add_segment_to_main_grid(main_grid, segment_grid)

    return main_grid

def _initialise_longitude_latitude_grid(
    spatial_bbox: tuple[float, float, float, float] = (-180.0, -90.0, 180.0, 90.0),
    spatial_grid_res: float = 0.05,
) -> xr.DataArray:
    r"""
    Create longitude-latitude grid of specified coordinates and spatial resolution.

    Parameters
    ----------
    spatial_bbox : tuple[float, float, float, float]
        Spatial bounding box, `(lon_min, lat_min, lon_max, lat_max)`, [:math:`\deg`]
    spatial_grid_res : float
        Spatial grid resolution, [:math:`\deg`]

    Returns
    -------
    xr.DataArray
        Longitude-latitude grid of specified coordinates and spatial resolution, filled with zeros.

    Notes
    -----
    This empty grid is used to store the aggregated plume properties of the individual
    contrail segments, such as the gridded contrail optical depth and radiative forcing.
    """
    lon_coords = np.arange(spatial_bbox[0], spatial_bbox[2] + spatial_grid_res, spatial_grid_res)
    lat_coords = np.arange(spatial_bbox[1], spatial_bbox[3] + spatial_grid_res, spatial_grid_res)
    return xr.DataArray(
        np.zeros((len(lon_coords), len(lat_coords))),
        dims=["longitude", "latitude"],
        coords={"longitude": lon_coords, "latitude": lat_coords},
    )

def segment_property_to_hi_res_grid(
    plume_segment: GeoVectorDataset,
    *,
    var_name: str,
    spatial_grid_res: float = 0.05,
    n_slices: int
) -> xr.DataArray:
    r"""
    Convert the plume segment property to a high-resolution longitude-latitude grid.

    Parameters
    ----------
    plume_segment : GeoVectorDataset
        Plume segment waypoints (head and tail).
    var_name : str
        Plume property of interest, where `var_name` must be included in `plume_segment`.
        For example, `tau_contrail`, `rf_sw`, `rf_lw`, and `rf_net`
    spatial_grid_res : float
        Spatial grid resolution, [:math:`\deg`]

    Returns
    -------
    xr.DataArray
        Plume segment dimension and property projected to a longitude-latitude grid.
    """
    # Ensure that `plume_segment` contains the required variables
    plume_segment.ensure_vars(("sin_a", "cos_a", "width", var_name))

    # Ensure that `plume_segment` only contains two waypoints and have the same time.
    assert len(plume_segment) == 2
    assert plume_segment["time"][0] == plume_segment["time"][1]

    # Calculate plume edges
    (
        plume_segment["lon_edge_l"],
        plume_segment["lat_edge_l"],
        plume_segment["lon_edge_r"],
        plume_segment["lat_edge_r"],
    ) = plume_edges(
        plume_segment["longitude"],
        plume_segment["latitude"],
        plume_segment["sin_a"],
        plume_segment["cos_a"],
        plume_segment["width"],
    )

    # Initialise plume segment grid with spatial domain that covers the plume area.
    lon_edges = np.concatenate(
        [plume_segment["lon_edge_l"], plume_segment["lon_edge_r"]], axis=0
    )
    lat_edges = np.concatenate(
        [plume_segment["lat_edge_l"], plume_segment["lat_edge_r"]], axis=0
    )

    spatial_bbox = geo.spatial_bounding_box(lon_edges, lat_edges, buffer=0.01)
    #print(spatial_bbox)
    segment_grid = _initialise_longitude_latitude_grid(spatial_bbox, spatial_grid_res)

    # # Get slice percentage and absolute mass in each segment and slice
    # slice_percentage = 1 / n_slices

    # segment_mass = ((plume_segment[var_name][0] 
    #                    + plume_segment[var_name][1]) / 2)
    

    # Calculate gridded plume segment properties
    for slice in range(n_slices):
        
        plume_slice = GeoVectorDataset()

        plume_slice["longitude"] = plume_segment["longitude"]
        plume_slice["latitude"] = plume_segment["latitude"]
        plume_slice["sin_a"] = plume_segment["sin_a"]
        plume_slice["cos_a"] = plume_segment["cos_a"]
        plume_slice["segment_width"] = plume_segment["width"]
        plume_slice["sigma_yy"] = plume_segment["sigma_yy"]
        plume_slice.attrs["slice_percentage"] = 1 / n_slices
        plume_slice.attrs["slice_mass"] = ((plume_segment[var_name][0] 
                                    + plume_segment[var_name][1]) / 2) / n_slices


        # Calculate plume slice lat lon positions
        lon_edges_slice, lat_edges_slice = plume_slices(plume_slice, slice)

        # Define shapely polygon
        slice_coords = ((lon_edges_slice[0], lat_edges_slice[0]),
                        (lon_edges_slice[1], lat_edges_slice[1]),
                        (lon_edges_slice[3], lat_edges_slice[3]),
                        (lon_edges_slice[2], lat_edges_slice[2]),
                        (lon_edges_slice[0], lat_edges_slice[0]))

        print(slice_coords)

        slice_polygon = Polygon(slice_coords)

        plt.plot(*slice_polygon.exterior.xy)
        plt.savefig("polygon")

        slice_area = slice_polygon.area

        print(slice_area)

        slice_grid = xr.DataArray(np.zeros_like(segment_grid), coords=segment_grid.coords, dims=segment_grid.dims)


        # Iterate over each cell in the grid
        cell_size = segment_grid.longitude[1] - segment_grid.longitude[0] # Grid cell size in degrees
        for i, lon in enumerate(segment_grid.longitude[:-1]):
            for j, lat in enumerate(segment_grid.latitude[:-1]):
                # Define the grid cell as a Shapely box
                cell = box(lon, lat, lon + cell_size, lat + cell_size)
                
                # Check intersection with the plume polygon
                if slice_polygon.intersects(cell):
                    intersection = slice_polygon.intersection(cell)
                    intersection_area = intersection.area
                    
                    # Store the intersection area in the grid
                    slice_grid[i, j] = (intersection_area / slice_area) * plume_slice.attrs["slice_mass"]

        segment_grid += slice_grid

        return(segment_grid)



def plume_edges(
    lon: npt.NDArray[np.float64],
    lat: npt.NDArray[np.float64],
    sin_a: npt.NDArray[np.float64],
    cos_a: npt.NDArray[np.float64],
    width: npt.NDArray[np.float64],
) -> tuple[
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
]:
    """
    Calculate the longitude and latitude of the plume edges to account for plume spreading.

    (lon_edge_l, lat_edge_l)        x---------------------

    (Plume midpoint: lon, lat)   X===================== ->

    (lon_edge_r, lat_edge_r)        x---------------------


    Parameters
    ----------
    lon : npt.NDArray[np.float64]
        longitude of plume waypoint, degrees
    lat : npt.NDArray[np.float64]
        latitude of plume waypoint, degrees
    sin_a : npt.NDArray[np.float64]
        sin(a), where a is the angle between the plume and the longitudinal axis
    cos_a : npt.NDArray[np.float64]
        cos(a), where a is the angle between the plume and the longitudinal axis
    width : npt.NDArray[np.float64]
        plume width at each waypoint, [:math:`m`]

    Returns
    -------
    tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]
        (lon_edge_l, lat_edge_l, lon_edge_r, lat_edge_r), longitudes and latitudes
        at the left and right edges of the plume, degrees
    """  # noqa: E501
    dlon = units.m_to_longitude_distance(width * sin_a * 0.5, lat)
    dlat = units.m_to_latitude_distance(width * cos_a * 0.5)

    lon_edge_l = lon - dlon
    lat_edge_l = lat + dlat
    lon_edge_r = lon + dlon
    lat_edge_r = lat - dlat

    return lon_edge_l, lat_edge_l, lon_edge_r, lat_edge_r

def plume_slices(
    plume_slice: GeoVectorDataset,
    slice: int,
) -> npt.NDArray[np.float64]:
    """
    Calculate the longitude and latitude of the plume slice points, to discretise the Gaussian distribution along the cross section of the plume.

    Parameters
    ----------

    Returns
    -------
    npt.NDArray[np.float64]
        (lon_slices, lat_slices), longitudes and latitudes
        at the increments of the plume cross section, degrees
    """  # noqa: E501
    std_dev = plume_slice["sigma_yy"] ** 0.5
    print(std_dev)

    lq = (plume_slice.attrs["slice_percentage"] / 4) + (plume_slice.attrs["slice_percentage"] / 2) * slice
    uq = 1 - lq

    mean = 0 # centreline of the plume

    z = norm.ppf(uq) - norm.ppf(lq)

    plume_slice["width"] = z*std_dev

    # Ensure slice_width does not exceed width
    plume_slice["width"] = np.minimum(plume_slice["width"], plume_slice["segment_width"])

    print(plume_slice["segment_width"])
    print(plume_slice["width"])

    # Calculate slice edges
    (
        plume_slice["lon_edge_l"],
        plume_slice["lat_edge_l"],
        plume_slice["lon_edge_r"],
        plume_slice["lat_edge_r"],
    ) = plume_edges(
        plume_slice["longitude"],
        plume_slice["latitude"],
        plume_slice["sin_a"],
        plume_slice["cos_a"],
        plume_slice["width"],
    )

    # Initialise plume segment grid with spatial domain that covers the plume area.
    lon_edges = np.concatenate(
        [plume_slice["lon_edge_l"], plume_slice["lon_edge_r"]], axis=0
    )
    lat_edges = np.concatenate(
        [plume_slice["lat_edge_l"], plume_slice["lat_edge_r"]], axis=0
    )

    return lon_edges, lat_edges

def _add_segment_to_main_grid(main_grid: xr.DataArray, segment_grid: xr.DataArray) -> xr.DataArray:
    """
    Add the gridded contrail segment to the main grid.

    Parameters
    ----------
    main_grid : xr.DataArray
        Aggregated contrail segment properties in a longitude-latitude grid.
    segment_grid : xr.DataArray
        Contrail segment dimension and property projected to a longitude-latitude grid.

    Returns
    -------
    xr.DataArray
        Aggregated contrail segment properties, including `segment_grid`.

    Notes
    -----
    - The spatial domain of `segment_grid` only covers the contrail segment, which is added to
        the `main_grid` which is expected to have a larger spatial domain than the `segment_grid`.
    - This architecture is used to reduce the computational resources.
    """
    lon_main = np.round(main_grid["longitude"].values, decimals=2)
    lat_main = np.round(main_grid["latitude"].values, decimals=2)

    lon_segment_grid = np.round(segment_grid["longitude"].values, decimals=2)
    lat_segment_grid = np.round(segment_grid["latitude"].values, decimals=2)

    main_grid_arr = main_grid.values
    subgrid_arr = segment_grid.values

    try:
        ix_ = np.searchsorted(lon_main, lon_segment_grid[0])
        ix = np.searchsorted(lon_main, lon_segment_grid[-1]) + 1
        iy_ = np.searchsorted(lat_main, lat_segment_grid[0])
        iy = np.searchsorted(lat_main, lat_segment_grid[-1]) + 1
    except IndexError:
        warnings.warn(
            "Contrail segment ignored as it is outside spatial bounding box of the main grid. "
        )
    else:
        main_grid_arr[ix_:ix, iy_:iy] = main_grid_arr[ix_:ix, iy_:iy] + subgrid_arr

    return xr.DataArray(main_grid_arr, coords=main_grid.coords)