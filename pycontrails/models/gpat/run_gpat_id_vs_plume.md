## Run GPAT for mass conservation validation test
So the base case scenario is the one to test all params against. So for each of the params to vary below, vary one at a time, while keeping everything else at base case.

Note that for hres_pl and hres_sim, these need to be automated at the same time, as I don't trust my interpolation scheme, and would rather keep them the same (but cba to recode it all).

also, i am running mine from sc_local.sh on my local box. Please see this and reconvert back to bc4 where necessary.

## Base case scenario
```python
# flight trajectory parameters
fl_params = {
    "t0_fl": pd.to_datetime("2022-01-20 13:00:00"),  # flight start time
    "rt_fl": pd.Timedelta(minutes=60),  # flight run time
    "ts_fl": pd.Timedelta(minutes=2),  # flight time step
    "ac_type": "A320",  # aircraft type
    "fl0_speed": 100.0,  # m/s
    "fl0_heading": 45.0,  # deg
    "fl0_coords0": (47.1, -32.9, 12500),  # lat, lon, alt [deg, deg, m]
    "sep_dist": (1000, 0, 0),  # dx, dy, dz [m]
    "n_ac": 2,  # number of aircraft
}

# plume dispersion parameters
plume_params = {
    "dt_integration": pd.Timedelta(minutes=2),  # integration time step
    "max_age": pd.Timedelta(hours=2),  # maximum age of the plume
    "depth": 50.0,  # initial plume depth, [m]
    "width": 50.0,  # initial plume width, [m]
    "shear": 0.01,  # wind shear [1/s]
    "hres_pl": 0.05, # horizontal resolution of the plume [deg]
    "vres_pl": 500, # vertical resolution of the plume [m]
    "n_slices": 10,  # number of plume slices
}

sim_params = {
    "t0_sim": pd.to_datetime("2022-01-20 12:00:00"),  # chemistry start time
    "rt_sim": plume_params["max_age"] + pd.Timedelta(hours=2),  # chemistry runtime
    "ts_sim": pd.Timedelta(seconds=20),  # chemistry time step
    "lat_bounds": (47.0, 48.0),  # lat bounds [deg]
    "lon_bounds": (-33.0, -32.0),  # lon bounds [deg]
    "alt_bounds": (12000, 13000),  # alt bounds [m]
    "hres_sim": 0.05,  # horizontal resolution [deg]
    "vres_sim": 500,  # vertical resolution [m]
    "eastward_wind": 0.0,  # m/s
    "northward_wind": 0.0,  # m/s
    "lagrangian_tendency_of_air_pressure": 0.0,  # m/s
    "species_in": np.array(["NO", "NO2", "CO", "HCHO", "CH3CHO", "C2H4", "C3H6", "C2H2", "BENZENE"]),
    "species_out": np.array(["O3", "NO2", "NO",
                            "NO3", "HNO3", "PAN",
                            "HONO", "HO2", "OH",
                            "H2O2", "CO", "HCHO",
                            "CH4"
                            ]),
    "job_id":   (f"mc_v_{fl_params['n_ac']}_"
                f"{fl_params['sep_dist'][0]}_"
                f"{fl_params['sep_dist'][1]}_"
                f"{plume_params['max_age'].components.hours}_"
                f"{plume_params["hres_pl"]}")
}
```
## Params to vary
- locations: NA, US, EU, SEA, SA

- fl_params["n_ac"]: [1, 2, 3, 5, 10]

- fl_params["sep_dist"][0]: [100, 1000, 2000, 5000, 10000] # dx [m]

- fl_params["sep_dist"][1]: [0, 100, 1000] # dy [m]

- plume_params["max_age"]: [1, 2, 5, 10, 12] # max age of plume waypoints [hours]

- plume_params["hres_pl"]: [0.01, 0.02, 0.05, 0.1, 0.5] # plume hres [degrees]
  plume_params["hres_sim"]: [0.01, 0.02, 0.05, 0.1, 0.5] # chem sim hres [degrees]
