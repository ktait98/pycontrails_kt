import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

pd.set_option("display.max_rows", 200)

from pycontrails import MetDataset

# from pycontrails.models.ps_model import PSFlight
# from pycontrails.models.emissions import Emissions
from pycontrails.ext.flight_gen import FlightGen
from pycontrails.models.boxmodel.boxm_ac import Boxm
from pycontrails.physics import units

def boxm_run(lons, lons_pl, lats, lats_pl, alts, times, fl_params, plume_params, chem_params, met_params):
        
        met = gen_met(lons, lats, alts, times, met_params)

        fl_gen = FlightGen(met, fl_params, plume_params, chem_params)

        fl = fl_gen.traj_gen()

        fl = fl_gen.calc_fb_emissions()

        fl_df, pl_df = fl_gen.sim_plumes()

        emi = fl_gen.plume_to_grid(lats_pl, lons_pl, alts, times)

        boxm = Boxm(met=met, fl_params=fl_params, plume_params=plume_params, chem_params=chem_params)
        
        boxm.eval(emi)

        return fl_df, pl_df, boxm

def gen_met(lons, lats, alts, times, met_params):

    # generate artifical met dataset (boxm currently only supports zero-wind scenarios)
    data_vars = {
        param: (
            ["longitude", "latitude", "level", "time"],
            da.full(
                (len(lons), len(lats), len(alts), len(times)),
                value,
                chunks=(len(lons), len(lats), len(alts), 100),
            ),
        )
        for param, value in met_params.items()
    }

    met = xr.Dataset(
        data_vars,
        coords={"longitude": lons, "latitude": lats, "level": units.m_to_pl(alts), "time": times},
    )

    met = MetDataset(met)

    return met
