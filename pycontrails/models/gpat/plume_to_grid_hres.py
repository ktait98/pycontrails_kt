"""Plume to grid module for aggregating plume segments to a high-resolution lat-lon grid."""

import warnings

import numpy as np
import numpy.typing as npt
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.ndimage import gaussian_filter
from shapely.geometry import Polygon, box

from pycontrails.core.vector import GeoVectorDataset
from pycontrails.physics import units
from pycontrails.utils import dependencies


def plume_to_grid_hres(
    time: pd.Timestamp | np.datetime64,
    plumes_t: GeoVectorDataset,
    *,
    var_name: str,
    main_grid: xr.DataArray,
    grid_res: tuple[float, float] = (0.05, 1000)

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
        "depth",
        "sigma_yy",
        "sigma_zz",
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

    # Add altitude and air pressure coordinates to the main grid
    main_grid = _add_vertical_coords(main_grid)

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
            pd.concat([heads_t[cols_req].loc[i], tails_t[cols_req].loc[i]], axis=1).T, copy=True)

        segment_grid = _segment_property_to_hi_res_grid(
            plume_segment, var_name=var_name, grid_res=grid_res)
        
        main_grid = _add_segment_to_main_grid(main_grid, segment_grid)

    return main_grid


def _initialise_3d_grid(
    spatial_bbox: tuple[float, float, float, float] = (-180.0, -90.0, 180.0, 90.0),
    grid_res: tuple[float, float] = (0.05, 1000),
) -> xr.DataArray:
    r"""
    Create 3-D grid of specified coordinates and spatial resolution.

    Parameters
    ----------
    spatial_bbox : tuple[float, float, float, float, float, float]
        Spatial bounding box, `(lon_min, lat_min, lon_max, lat_max, alt_min, alt_max)`, [:math:`\deg`, :math:`m`]
    grid_res : tuple[float, float]
        Grid resolution, `(horiz_res, vert_res) [:math:`\deg`, :math:`m`]`

    Returns
    -------
    xr.DataArray
        Longitude-latitude-altitude grid of specified coordinates and spatial resolution, filled with zeros.

    Notes
    -----
    This empty grid is used to store the aggregated plume properties of the individual
    contrail segments, such as the gridded plume optical depth and radiative forcing.
    """
    lon_coords = np.arange(spatial_bbox[0], spatial_bbox[2] + grid_res[0], grid_res[0])
    lat_coords = np.arange(spatial_bbox[1], spatial_bbox[3] + grid_res[0], grid_res[0])
    alt_coords = np.arange(spatial_bbox[4], spatial_bbox[5] + grid_res[1], grid_res[1])

    return xr.DataArray(
        np.zeros((len(lon_coords), len(lat_coords), len(alt_coords))),
        dims=["longitude", "latitude", "altitude"],
        coords={"longitude": lon_coords, "latitude": lat_coords, "altitude": alt_coords}
    )


def _segment_property_to_hi_res_grid(
    plume_segment: GeoVectorDataset, *, var_name: str,
    grid_res: tuple[float, float] = (0.05, 1000)
) -> xr.DataArray:
    r"""
    Convert the plume segment property to a high-resolution longitude-latitude-altitude grid.

    Parameters
    ----------
    plume_segment : GeoVectorDataset
        Plume segment waypoints (head and tail).
    var_name : str
        Plume property of interest, where `var_name` must be included in `plume_segment`.
        For example, `tau_contrail`, `rf_sw`, `rf_lw`, and `rf_net`
    coarse_grid_res : tuple[float, float]
        Coarse grid resolution, `(coarse_hres, coarse_vres)`, [:math:`deg`, :math:`m`]
    fine_grid_res : tuple[float, float]
        Fine grid resolution, `(fine_hres, fine_vres)`, [:math:`\deg`, :math:`m`]

    Returns
    -------
    xr.DataArray
        Plume segment dimension and property projected to a 3D grid.
    """
    # Ensure that `plume_segment` contains the required variables
    plume_segment.ensure_vars(("sin_a", "cos_a", "width", "depth", "sigma_yy", "sigma_zz", var_name))

    # Ensure that `plume_segment` only contains two waypoints and have the same time.
    assert len(plume_segment) == 2
    assert plume_segment["time"][0] == plume_segment["time"][1]

    # Calculate plume edges
    (
        plume_segment["lon_edge_l"],
        plume_segment["lat_edge_l"],
        plume_segment["alt_edge_u"],
        plume_segment["lon_edge_r"],
        plume_segment["lat_edge_r"],
        plume_segment["alt_edge_d"]
    ) = plume_edges(
        plume_segment["longitude"],
        plume_segment["latitude"],
        plume_segment["altitude"],
        plume_segment["sin_a"],
        plume_segment["cos_a"],
        plume_segment["width"],
        plume_segment["depth"],
    )

    # Initialise plume segment grid with spatial domain that covers the plume area.
    lon_edges = np.concatenate([plume_segment["lon_edge_l"], plume_segment["lon_edge_r"]], axis=0)
    lat_edges = np.concatenate([plume_segment["lat_edge_l"], plume_segment["lat_edge_r"]], axis=0)
    alt_edges = np.concatenate([plume_segment["alt_edge_u"], plume_segment["alt_edge_d"]], axis=0)

    # Create a spatial bounding box for the plume segment
    fine_spatial_bbox = _spatial_bounding_box(
        lon_edges, lat_edges, alt_edges, grid_res, buffer=grid_res
    )
    
    # Create a 3D grid of specified coordinates and spatial resolution
    segment_grid = _initialise_3d_grid(fine_spatial_bbox, grid_res)

    # Calculate weights for the segment grid
    weights = _pixel_weights(plume_segment, segment_grid)

    # Calculate the shortest distance from the plume segment to each pixel in the grid
    dist_perp = _segment_perpendicular_distance_to_pixels(plume_segment, weights)

    # Calculate the plume concentration at each pixel
    plume_concentration = _gaussian_plume_concentration(
        plume_segment,
        weights,
        dist_perp,
    )
    
    # Distribute selected contrail property to grid
    plume_property = plume_concentration * (
        weights * xr.ones_like(weights) * plume_segment[var_name][1]
        + (1 - weights) * xr.ones_like(weights) * plume_segment[var_name][0]
    )

    # Smooth the plume concentration using a Gaussian filter
    plume_property = _gaussian_filter_plume_property(plume_property, sigma=0.5)
    
    # Normalise the plume concentration to ensure mass conservation
    plume_property = _normalise_plume_property(plume_property, plume_segment)
   
    return plume_property


def plume_edges(
    lon: npt.NDArray[np.float64],
    lat: npt.NDArray[np.float64],
    alt: npt.NDArray[np.float64],
    sin_a: npt.NDArray[np.float64],
    cos_a: npt.NDArray[np.float64],
    width: npt.NDArray[np.float64],
    depth: npt.NDArray[np.float64],
) -> tuple[
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
]:
    """Calculate the longitude, latitude and altitude of the plume edges to account for plume spreading.

    (lon_edge_l, lat_edge_l)        x---------------------

    (Plume midpoint: lon, lat)   X===================== ->

    (lon_edge_r, lat_edge_r)        x---------------------

    Parameters
    ----------
    lon : npt.NDArray[np.float64]
        longitude of plume waypoint, degrees
    lat : npt.NDArray[np.float64]
        latitude of plume waypoint, degrees
    alt : npt.NDArray[np.float64]
        altitude of plume waypoint, [:math:`m`]
    sin_a : npt.NDArray[np.float64]
        sin(a), where a is the angle between the plume and the longitudinal axis
    cos_a : npt.NDArray[np.float64]
        cos(a), where a is the angle between the plume and the longitudinal axis
    width : npt.NDArray[np.float64]
        plume width at each waypoint, [:math:`m`]
    depth : npt.NDArray[np.float64]
        plume depth at each waypoint, [:math:`m`]

    Returns
    -------
    tuple[npt.NDArray[np.float64], 
          npt.NDArray[np.float64], 
          npt.NDArray[np.float64], 
          npt.NDArray[np.float64], 
          npt.NDArray[np.float64], 
          npt.NDArray[np.float64]]
        (lon_edge_l, lat_edge_l, lon_edge_r, lat_edge_r, alt_edge_u, alt_edge_d), 
        lons, lats and alts at the edges of the plume [degrees, degrees, m]
    """
    dlon = units.m_to_longitude_distance(width * sin_a * 0.5, lat)
    dlat = units.m_to_latitude_distance(width * cos_a * 0.5)
    dalt = depth * 0.5

    lon_edge_l = lon - dlon
    lat_edge_l = lat + dlat
    lon_edge_r = lon + dlon
    lat_edge_r = lat - dlat
    alt_edge_u = alt + dalt
    alt_edge_d = alt - dalt

    return lon_edge_l, lat_edge_l, alt_edge_u, lon_edge_r, lat_edge_r, alt_edge_d


def _pixel_weights(plume_segment: GeoVectorDataset, segment_grid: xr.DataArray) -> xr.DataArray:
    head = plume_segment.dataframe.iloc[0]
    tail = plume_segment.dataframe.iloc[1]

    # Calculate determinant
    dx = units.longitude_distance_to_m(
        (tail["longitude"] - head["longitude"]),
        0.5 * (head["latitude"] + tail["latitude"]),
    )
    dy = units.latitude_distance_to_m(tail["latitude"] - head["latitude"])
    dz = tail["altitude"] - head["altitude"]
    det = dx**2 + dy**2 + dz**2
    if det == 0:
        raise ValueError("Plume segment has zero length.")
    if det < 0:
        raise ValueError("Plume segment has negative length.")

    # Use indexing='ij' for correct axis order
    lon_grid, lat_grid, alt_grid = np.meshgrid(
        segment_grid["longitude"].values,
        segment_grid["latitude"].values,
        segment_grid["altitude"].values,
        indexing='ij'
    )

    dx_grid = units.longitude_distance_to_m(
        (lon_grid - head["longitude"]),
        0.5 * (head["latitude"] + lat_grid),
    )
    dy_grid = units.latitude_distance_to_m(lat_grid - head["latitude"])
    dz_grid = alt_grid - head["altitude"]

    weights = (dx * dx_grid + dy * dy_grid + dz * dz_grid) / det

    return xr.DataArray(
        data=weights,
        dims=["longitude", "latitude", "altitude"],
        coords={
            "longitude": segment_grid["longitude"],
            "latitude": segment_grid["latitude"],
            "altitude": segment_grid["altitude"],
        },
    )

def _segment_perpendicular_distance_to_pixels(
    plume_segment: GeoVectorDataset, weights: xr.DataArray
) -> xr.DataArray:
    head = plume_segment.dataframe.iloc[0]
    tail = plume_segment.dataframe.iloc[1]

    # Use indexing='ij' for correct axis order
    lon_grid, lat_grid, alt_grid = np.meshgrid(
        weights["longitude"].values,
        weights["latitude"].values,
        weights["altitude"].values,
        indexing='ij'
    )

    lon_s = head["longitude"] + weights.values * (tail["longitude"] - head["longitude"])
    lat_s = head["latitude"] + weights.values * (tail["latitude"] - head["latitude"])
    alt_s = head["altitude"] + weights.values * (tail["altitude"] - head["altitude"])

    lon_dist = units.longitude_distance_to_m(np.abs(lon_grid - lon_s), 0.5 * (lat_s + lat_grid))
    lat_dist = units.latitude_distance_to_m(np.abs(lat_grid - lat_s))
    alt_dist = np.abs(alt_grid - alt_s)

    dist_perp_h = (lon_dist**2 + lat_dist**2) ** 0.5
    dist_perp_v = np.abs(alt_dist)

    dist_perp = np.stack([dist_perp_h, dist_perp_v], axis=0)
    
    return xr.DataArray(
    dist_perp,
    dims=["perp_type", "longitude", "latitude", "altitude"],
    coords={
        "perp_type": ["horizontal", "vertical"],
        "longitude": weights["longitude"],
        "latitude": weights["latitude"],
        "altitude": weights["altitude"],
    },
)


def _gaussian_plume_concentration(
    plume_segment: GeoVectorDataset,
    weights: xr.DataArray,
    dist_perpendicular: xr.DataArray,
) -> xr.DataArray:
    """
    Calculate relative gaussian plume concentration along the plume width.

    Parameters
    ----------
    plume_segment : GeoVectorDataset
        Plume segment waypoints (head and tail).
    weights : xr.DataArray
        Pixel weights for `segment_grid`.
        See `_pixel_weights` function.
    dist_perpendicular : xr.DataArray
        Perpendicular distance from plume segment to each segment grid pixel, [:math:`m`]
        See `_segment_perpendicular_distance_to_pixels` function.

    Returns
    -------
    xr.DataArray
        Relative gaussian plume concentration along the plume width

    Notes
    -----
    - 2D Gaussian plume concentration is calculated using the following equation:
        .. math::

    - See Appendix A11 of :cite:`schumannContrailCirrusPrediction2012`.
    """
    head = plume_segment.dataframe.iloc[0]
    tail = plume_segment.dataframe.iloc[1]

    sigma_yy = weights.values * tail["sigma_yy"] + (1 - weights.values) * head["sigma_yy"]
    sigma_zz = weights.values * tail["sigma_zz"] + (1 - weights.values) * head["sigma_zz"]
    
    dist_perp_h = dist_perpendicular.sel(perp_type="horizontal").values
    dist_perp_v = dist_perpendicular.sel(perp_type="vertical").values

    concentration = np.where(
        (weights.values < 0) | (weights.values > 1),
        0,
        (1 / (2 * np.pi * sigma_yy * sigma_zz) ** 0.5) *  # Normalization factor
        np.exp(
            -0.5 * (
                (dist_perp_h**2 / sigma_yy) +  # Horizontal distance term
                (dist_perp_v**2 / sigma_zz)   # Vertical distance term
            )
        )
    )

    return xr.DataArray(concentration, coords=weights.coords)


def _gaussian_filter_plume_property(plume_property: xr.DataArray, sigma: float) -> xr.DataArray:
    """
    Smooth the property grid using a Gaussian filter.

    Parameters
    ----------
    plume_property : xr.DataArray
        The input 3D grid (longitude, latitude, altitude).
    sigma : float
        The standard deviation of the Gaussian kernel.

    Returns
    -------
    xr.DataArray
        The smoothed grid.
    """
    # Apply Gaussian filter
    smoothed_values = gaussian_filter(plume_property.values, sigma=sigma, mode="constant", cval=0.0)

    # Return the smoothed grid as an xarray DataArray
    return xr.DataArray(
        smoothed_values,
        dims=plume_property.dims,
        coords=plume_property.coords,
    )


def _normalise_plume_property(
    plume_property: xr.DataArray,
    plume_segment: GeoVectorDataset,
) -> xr.DataArray:
    """
    Normalise the plume segment grid to ensure mass conservation.

    Parameters
    ----------
    plume_property : xr.DataArray
        Plume property grid with spatial domain that covers the plume area.
    plume_segment : GeoVectorDataset
        Plume segment waypoints (head and tail).

    Returns
    -------
    xr.DataArray
        Normalised plume segment grid.
    """
    # Calculate the total mass of the plume segment
    total_mass = np.sum(plume_segment["width"] * plume_segment["depth"])

    # Calculate the total mass in the segment grid
    total_mass_grid = np.sum(plume_property.values)

    # Normalise the segment grid to ensure mass conservation
    if total_mass_grid > 0:
        normalised_plume_property = (plume_property / total_mass_grid) * total_mass
    else:
        normalised_plume_property = plume_property

    return normalised_plume_property


def _add_segment_to_main_grid(main_grid: xr.DataArray, segment_grid: xr.DataArray) -> xr.DataArray:
    r"""
    Add the gridded plume segment to the main grid.

    Parameters
    ----------
    main_grid : xr.DataArray
        Aggregated plume segment properties in a longitude-latitude grid.
    segment_grid : xr.DataArray
        plume segment dimension and property projected to a longitude-latitude grid.

    Returns
    -------
    xr.DataArray
        Aggregated plume segment properties, including `segment_grid`.

    Notes
    -----
    - The spatial domain of `segment_grid` only covers the plume segment, which is added to
        the `main_grid` which is expected to have a larger spatial domain than the `segment_grid`.
    - This architecture is used to reduce the computational resources.
    """


    lon_main = np.round(main_grid["longitude"].values, decimals=2)
    lat_main = np.round(main_grid["latitude"].values, decimals=2)
    alt_main = np.round(main_grid["altitude"].values, decimals=2)

    lon_segment_grid = np.round(segment_grid["longitude"].values, decimals=2)
    lat_segment_grid = np.round(segment_grid["latitude"].values, decimals=2)
    alt_segment_grid = np.round(segment_grid["altitude"].values, decimals=2)

    main_grid_arr = main_grid.values
    subgrid_arr = segment_grid.values

    try:
        ix_ = np.searchsorted(lon_main, lon_segment_grid[0])
        ix = np.searchsorted(lon_main, lon_segment_grid[-1]) + 1
        iy_ = np.searchsorted(lat_main, lat_segment_grid[0])
        iy = np.searchsorted(lat_main, lat_segment_grid[-1]) + 1
        iz_ = np.searchsorted(alt_main, alt_segment_grid[0])
        iz = np.searchsorted(alt_main, alt_segment_grid[-1]) + 1

        # Create a mask to ensure the subgrid fits within the main grid
        mask = np.zeros_like(main_grid_arr[ix_:ix, iy_:iy, iz_:iz])
        mask[: subgrid_arr.shape[0], : subgrid_arr.shape[1], : subgrid_arr.shape[2]] = subgrid_arr

        # Add the masked subgrid to the main grid
        main_grid_arr[ix_:ix, iy_:iy, iz_:iz] += mask

    except (IndexError, ValueError) as e:
        warnings.warn(f"plume segment resized due to {e}. ")

    return xr.DataArray(main_grid_arr, coords=main_grid.coords)


def _spatial_bounding_box(
    longitude: npt.NDArray[np.float64],
    latitude: npt.NDArray[np.float64],
    altitude: npt.NDArray[np.float64],
    spatial_grid_res: tuple[float, float] = (0.5, 1000.0),
    buffer: tuple[float, float] = (0.1, 100.0),
) -> tuple[float, float, float, float]:
    r"""
    Construct rectangular spatial bounding box from a set of waypoints.

    Parameters
    ----------
    longitude : np.ndarray
        1D Longitude values with index corresponding to longitude inputs, [:math:`\deg`]
    latitude : np.ndarray
        1D Latitude values with index corresponding to latitude inputs, [:math:`\deg`]
    spatial_grid_res: float
        Horiz grid res that rounds the corner positions to the nearest grid cell edge,
        [:math:`\deg`]
    buffer: float
        Add buffer to rectangular spatial bounding box, [:math:`\deg`]

    Returns
    -------
    tuple[float, float, float, float]
        Spatial bounding box, ``(lon_min, lat_min, lon_max, lat_max)``, [:math:`\deg`]

    Examples
    --------
    >>> rng = np.random.default_rng(654321)
    >>> lon = rng.uniform(-180, 180, size=30)
    >>> lat = rng.uniform(-90, 90, size=30)
    >>> spatial_bounding_box(lon, lat)
    (np.float64(-168.0), np.float64(-77.0), np.float64(155.0), np.float64(82.0))
    """
    lon_min = max((np.min(longitude) - buffer[0]), -180.0)
    lon_max = min((np.max(longitude) + buffer[0]), 179.99)
    lat_min = max((np.min(latitude) - buffer[0]), -90.0)
    lat_max = min((np.max(latitude) + buffer[0]), 90.0)
    alt_min = max((np.min(altitude) - buffer[1]), 0.0)
    alt_max = min((np.max(altitude) + buffer[1]), 60000.0)

    lon_min = round(lon_min / spatial_grid_res[0]) * spatial_grid_res[0]
    lon_max = round(lon_max / spatial_grid_res[0]) * spatial_grid_res[0]
    lat_min = round(lat_min / spatial_grid_res[0]) * spatial_grid_res[0]
    lat_max = round(lat_max / spatial_grid_res[0]) * spatial_grid_res[0]
    alt_min = round(alt_min / spatial_grid_res[1]) * spatial_grid_res[1]
    alt_max = round(alt_max / spatial_grid_res[1]) * spatial_grid_res[1]
    return lon_min, lat_min, lon_max, lat_max, alt_min, alt_max


def _add_vertical_coords(data: xr.Dataset) -> xr.Dataset:
    """Add "air_pressure" and "altitude" coordinates to data.

    .. versionchanged:: 0.52.1
        Ensure that the ``dtype`` of the additional vertical coordinates agree
        with the ``dtype`` of the underlying gridded data.
    """
    data["level"].attrs.update(units="hPa", long_name="Pressure", positive="down")

    # XXX: use the dtype of the data to determine the precision of these coordinates
    # There are two competing conventions here:
    # - coordinate data should be float64
    # - gridded data is typically float32
    # - air_pressure and altitude often play both roles
    # It is more important for air_pressure and altitude to be grid-aligned than to be
    # coordinate-aligned, so we use the dtype of the data to determine the precision of
    # these coordinates
    level = data["level"].values

    if "air_pressure" not in data.coords:
        data = data.assign_coords(air_pressure=("level", level * 100.0))

    if "altitude" not in data.coords:
        data = data.assign_coords(altitude=("level", units.pl_to_m(level)))

    return data