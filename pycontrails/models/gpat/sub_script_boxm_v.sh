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

# ./run_gpat_boxm_v.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --run_gpat true --job_id "boxm_v_NA_20_5"
# ./run_gpat_boxm_v.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "6000,7000" --t0_sim "2022-11-10T12:00:00" --run_gpat true --job_id "boxm_v_US_20_5"
./run_gpat_boxm_v.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "8000,9000" --t0_sim "2022-05-20T12:00:00" --run_gpat true --job_id "boxm_v_EU_20_5"
./run_gpat_boxm_v.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "2000,3000" --t0_sim "2022-03-05T12:00:00" --run_gpat true --job_id "boxm_v_SEA_20_5"
./run_gpat_boxm_v.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "2022-08-15T12:00:00" --run_gpat true --job_id "boxm_v_SA_20_5"

#./run_gpat_boxm_v.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --run_gpat true --job_id "boxm_v_NA_20_5"

