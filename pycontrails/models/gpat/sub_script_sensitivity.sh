#!/bin/bash

#SBATCH --job-name=BoxModels
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:0:0
#SBATCH --mem=64G
#SBATCH --account=aero004481

#SBATCH --ntasks-per-node=1

#Apparently this requires the entire memory on the node, which is what we want rather than --exclusive

#Don't care if others use the node but use no memory on it! Tricky to manage probably.

#SBATCH --mem=0

#aero004301 
#aero004481 
   #default 
#isys015562

conda activate contrails

NCLIBS=`nc-config --libdir`
NFLIBS=`nf-config --prefix`/lib


export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NCLIBS}:${NFLIBS}
export PYCONTRAILSDIR=/user/work/${USER}/pycontrails_kt/pycontrails/

# ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --n_ac 1 --max_age 2 --run_gpat --job_id "sensitivity_NA_1_1000_0_2_0.01_0.05"
# ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --n_ac 2 --max_age 2 --run_gpat --job_id "sensitivity_NA_2_1000_0_2_0.01_0.05"
# ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --n_ac 5 --max_age 2 --run_gpat --job_id "sensitivity_NA_5_1000_0_2_0.01_0.05"

./run_gpat_sensitivity.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "11000,12000" --t0_sim "2022-11-10T12:00:00" \
                          --fl0_coords0 "37.1,-96.9,11500" --t0_fl "2022-11-10T13:00:00" --n_ac 1 --max_age 2 --run_gpat --job_id "sensitivity_US_1_1000_0_2_0.01_0.05"
./run_gpat_sensitivity.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "11000,12000" --t0_sim "2022-11-10T12:00:00" \
                          --fl0_coords0 "37.1,-96.9,11500" --t0_fl "2022-11-10T13:00:00" --n_ac 2 --max_age 2 --run_gpat --job_id "sensitivity_US_2_1000_0_2_0.01_0.05"
./run_gpat_sensitivity.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "11000,12000" --t0_sim "2022-11-10T12:00:00" \
                          --fl0_coords0 "37.1,-96.9,11500" --t0_fl "2022-11-10T13:00:00" --n_ac 5 --max_age 2 --run_gpat --job_id "sensitivity_US_5_1000_0_2_0.01_0.05"

./run_gpat_sensitivity.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "9000,10000" --t0_sim "2022-05-20T12:00:00" \
                          --fl0_coords0 "42.1,7.1,9500" --t0_fl "2022-05-20T13:00:00" --n_ac 1 --max_age 2 --run_gpat --job_id "sensitivity_EU_1_1000_0_2_0.01_0.05"
./run_gpat_sensitivity.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "9000,10000" --t0_sim "2022-05-20T12:00:00" \
                          --fl0_coords0 "42.1,7.1,9500" --t0_fl "2022-05-20T13:00:00" --n_ac 2 --max_age 2 --run_gpat --job_id "sensitivity_EU_2_1000_0_2_0.01_0.05"
./run_gpat_sensitivity.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "9000,10000" --t0_sim "2022-05-20T12:00:00" \
                          --fl0_coords0 "42.1,7.1,9500" --t0_fl "2022-05-20T13:00:00" --n_ac 5 --max_age 2 --run_gpat --job_id "sensitivity_EU_5_1000_0_2_0.01_0.05"

./run_gpat_sensitivity.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim "2022-03-05T12:00:00" \
                          --fl0_coords0 "22.1,102.1,10500" --t0_fl "2022-03-05T13:00:00" --n_ac 1 --max_age 2 --run_gpat --job_id "sensitivity_SEA_1_1000_0_2_0.01_0.05"
./run_gpat_sensitivity.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim "2022-03-05T12:00:00" \
                          --fl0_coords0 "22.1,102.1,10500" --t0_fl "2022-03-05T13:00:00" --n_ac 2 --max_age 2 --run_gpat --job_id "sensitivity_SEA_2_1000_0_2_0.01_0.05"
./run_gpat_sensitivity.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim "2022-03-05T12:00:00" \
                          --fl0_coords0 "22.1,102.1,10500" --t0_fl "2022-03-05T13:00:00" --n_ac 5 --max_age 2 --run_gpat --job_id "sensitivity_SEA_5_1000_0_2_0.01_0.05"

./run_gpat_sensitivity.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "2022-08-15T12:00:00" \
                          --fl0_coords0 " -27.9,-67.9,13500" --t0_fl "2022-08-15T13:00:00" --n_ac 1 --max_age 2 --run_gpat --job_id "sensitivity_SA_1_1000_0_2_0.01_0.05"
./run_gpat_sensitivity.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "2022-08-15T12:00:00" \
                          --fl0_coords0 " -27.9,-67.9,13500" --t0_fl "2022-08-15T13:00:00" --n_ac 2 --max_age 2 --run_gpat --job_id "sensitivity_SA_2_1000_0_2_0.01_0.05"
./run_gpat_sensitivity.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "2022-08-15T12:00:00" \
                          --fl0_coords0 " -27.9,-67.9,13500" --t0_fl "2022-08-15T13:00:00" --n_ac 5 --max_age 2 --run_gpat --job_id "sensitivity_SA_5_1000_0_2_0.01_0.05"