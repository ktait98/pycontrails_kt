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

# bg_run_<loc>_<month>

# datetimes=("2022-05-01T12:00:00" "2022-06-01T12:00:00" "2022-07-01T12:00:00")

# ./run_gpat_bg_runs.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim "2022-04-01T12:00:00" --run_gpat --job_id "bg_run_SEA_2022_04_01T12_00_00"
# ./run_gpat_bg_runs.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "2022-04-01T12:00:00" --run_gpat --job_id "bg_run_SA_2022_04_01T12_00_00"

# for datetime in "${datetimes[@]}"; do
#    # Replace colons with underscores
#    datetime_safe=$(echo $datetime | tr ':' '_')
   
#    ./run_gpat_bg_runs.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim $datetime --run_gpat --job_id "bg_run_NA_${datetime_safe}"
#    ./run_gpat_bg_runs.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "11000,12000" --t0_sim $datetime --run_gpat --job_id "bg_run_US_${datetime_safe}"
#    ./run_gpat_bg_runs.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "9000,10000" --t0_sim $datetime --run_gpat --job_id "bg_run_EU_${datetime_safe}"
#    ./run_gpat_bg_runs.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim $datetime --run_gpat --job_id "bg_run_SEA_${datetime_safe}"
#    ./run_gpat_bg_runs.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim $datetime --run_gpat --job_id "bg_run_SA_${datetime_safe}"

# done

# ./run_gpat_bg_runs.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "11000,12000" --t0_sim "2022-07-01T12:00:00" --run_gpat --job_id "bg_run_US_2022_07_01T12_00_00"
# ./run_gpat_bg_runs.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "11000,12000" --t0_sim "2022-10-01T12:00:00" --run_gpat --job_id "bg_run_US_2022_10_01T12_00_00"


# ./run_gpat_bg_runs.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "9000,10000" --t0_sim "2022-07-01T12:00:00" --run_gpat --job_id "bg_run_EU_07_01T12_00_00"
# ./run_gpat_bg_runs.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "9000,10000" --t0_sim "2022-10-01T12:00:00" --run_gpat --job_id "bg_run_EU_10_01T12_00_00"


# ./run_gpat_bg_runs.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-07-01T12:00:00" --run_gpat --job_id "bg_run_NA_07_01T12_00_00"
# ./run_gpat_bg_runs.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-10-01T12:00:00" --run_gpat --job_id "bg_run_NA_10_01T12_00_00"


# ./run_gpat_bg_runs.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim "2022-07-01T12:00:00" --run_gpat --job_id "bg_run_SEA_07_01T12_00_00"


# ./run_gpat_bg_runs.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "2022-07-01T12:00:00" --run_gpat --job_id "bg_run_SA_07_01T12_00_00"
# ./run_gpat_bg_runs.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "2022-10-01T12:00:00" --run_gpat --job_id "bg_run_SA_10_01T12_00_00"

./run_gpat_bg_runs.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim "2022-10-01T12:00:00" --run_gpat --job_id "bg_run_SEA_07_01T12_00_00"