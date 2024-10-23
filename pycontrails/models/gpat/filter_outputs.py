import os
import numpy as np
import pandas as pd
import pickle

def load_params(job_id):
    params_path = f"outputs/{job_id}/params_{job_id}.pkl"
    with open(params_path, 'rb') as f:
        params = pd.read_pickle("outputs/" + job_id + "/params_" + job_id + ".pkl")
        return params
    
def match_params(params, criteria):
    for outer_key, inner_criteria in criteria.items():
        if outer_key not in params:
            return False
        for key, value in inner_criteria.items():
            if key not in params[outer_key] or params[outer_key][key] != value:
                return False
    return True

def search_jobs(criteria):
    matching_job_ids = []
    outputs_dir = "outputs"

    for job_id in os.listdir(outputs_dir):
        job_dir = os.path.join(outputs_dir, job_id)
        if os.path.isdir(job_dir):
            try:
                params = load_params(job_id)
                if match_params(params, criteria):
                    matching_job_ids.append(job_id)
            except Exception as e:
                print(f"Error loading params for job {job_id}: {e}")

    return matching_job_ids

# Example criteria
criteria = {
    "sim_params": {
        "rt_sim": pd.Timedelta(hours=12),  # chemistry runtime
    },
    "fl_params": {
        "n_ac": 1,  # number of aircraft
    }
}

matching_job_ids = search_jobs(criteria)
print("Matching job IDs:", matching_job_ids)