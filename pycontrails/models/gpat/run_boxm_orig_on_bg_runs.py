#!/usr/bin/env python
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
import os
import re
from pycontrails.core import GeoVectorDataset
from pycontrails.models.gpat.gpat import mc_test, boxm_test
from pycontrails.models.gpat.pp_gpat import GPATPostProcessor

outputs_dir = "/projects/Impact_of_aviation_on_climate/Kieran2024/outputs/" 

# Filter criteria
criteria = {"job_id": "bg_run"}

pp_gpat = GPATPostProcessor(outputs_dir, criteria)

# Create dicts to hold all necessary data (but no more)
chem_ds_dict = {}

# Check if the dictionary file exists
if os.path.exists(outputs_dir + 'chem_ds_dict.pkl'):
    # Load the dictionary from the file
    with open(outputs_dir + 'chem_ds_dict.pkl', 'rb') as f:
        chem_ds_dict = pickle.load(f)
    print("Loaded chem_ds_dict from file.")
else:
    # Initialize the dictionary
    chem_ds_dict = {}
    
    # Load the data
    for job_id in pp_gpat.job_ids:
        chem_ds_dict[job_id] = pp_gpat.load_chem_ds(job_id)

    # Save the dictionary to a file
    with open(outputs_dir + 'chem_ds_dict.pkl', 'wb') as f:
        pickle.dump(chem_ds_dict, f)
    print("Saved chem_ds_dict to file.")

run_path = "/user/work/kt16229/pycontrails_kt/pycontrails/models/gpat/"
data_path = "/projects/Impact_of_aviation_on_climate/Kieran2024/"
# Loop through each job and collect the data
for job_id, chem_ds in chem_ds_dict.items():

    # Define chem ds as central cell (point analysis)
    chem_ds = pp_gpat.calc_cell(chem_ds, 1, 1)
    chem_ds = chem_ds.isel(job_id=0)

    chem_ds_dict[job_id] = boxm_test(run_path, data_path, job_id, chem_ds)
    print("Finished boxm_test for " + job_id)