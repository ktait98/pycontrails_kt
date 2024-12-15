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

chmod -R 755 /projects/Impact_of_aviation_on_climate/Kieran2024/

loc=$1
if [ $loc == "NA" ]; then
   datetimes=("2022-07-01T12:00:00" "2022-12-01T12:00:00")
   t0_fl=("2022-07-01T13:00:00" "2022-12-01T13:00:00")

   for i in "${!datetimes[@]}"; do
      datetime=${datetimes[$i]}
      t0_fl_value=${t0_fl[$i]}

      # Replace colons with underscores
      datetime_safe=$(echo $datetime | tr ':' '_')

      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                                --n_ac 1 --fl0_coords0 "47.5,-32.9,12500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_NA_1_${datetime_safe}"
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                                --n_ac 5 --fl0_coords0 "47.5,-32.9,12500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_NA_5_${datetime_safe}"
   done

elif [ $loc == "US" ]; then
   datetimes=("2022-05-01T12:00:00" "2022-01-01T12:00:00")
   t0_fl=("2022-05-01T13:00:00" "2022-01-01T13:00:00")

   for i in "${!datetimes[@]}"; do
      datetime=${datetimes[$i]}
      t0_fl_value=${t0_fl[$i]}

      # Replace colons with underscores
      datetime_safe=$(echo $datetime | tr ':' '_')

      ./run_gpat_sensitivity.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "11000,12000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                              --n_ac 1 --fl0_coords0 "37.5,-96.9,11500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_US_1_${datetime_safe}"
      ./run_gpat_sensitivity.py --lat_bounds "37.0,38.0" --lon_bounds " -97.0,-96.0" --alt_bounds "11000,12000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                              --n_ac 5 --fl0_coords0 "37.5,-96.9,11500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_US_5_${datetime_safe}"
   done

elif [ $loc == "EU" ]; then
   datetimes=("2022-06-01T12:00:00" "2022-12-01T12:00:00")
   t0_fl=("2022-06-01T13:00:00" "2022-12-01T13:00:00")
      
   for i in "${!datetimes[@]}"; do
      datetime=${datetimes[$i]}
      t0_fl_value=${t0_fl[$i]}

      # Replace colons with underscores
      datetime_safe=$(echo $datetime | tr ':' '_')

      ./run_gpat_sensitivity.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "9000,10000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                              --n_ac 1 --fl0_coords0 "42.5,7.1,9500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_EU_1_${datetime_safe}"
      ./run_gpat_sensitivity.py --lat_bounds "42.0,43.0" --lon_bounds "7.0,8.0" --alt_bounds "9000,10000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                              --n_ac 5 --fl0_coords0 "42.5,7.1,9500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_EU_5_${datetime_safe}"
   done

elif [ $loc == "SEA" ]; then
   datetimes=("2022-02-01T12:00:00" "2022-06-01T12:00:00")
   t0_fl=("2022-02-01T13:00:00" "2022-06-01T13:00:00")
      
   for i in "${!datetimes[@]}"; do
      datetime=${datetimes[$i]}
      t0_fl_value=${t0_fl[$i]}

      # Replace colons with underscores
      datetime_safe=$(echo $datetime | tr ':' '_')

      ./run_gpat_sensitivity.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                              --n_ac 1 --fl0_coords0 "22.5,102.1,10500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_SEA_1_${datetime_safe}"
      ./run_gpat_sensitivity.py --lat_bounds "22.0,23.0" --lon_bounds "102.0,103.0" --alt_bounds "10000,11000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                              --n_ac 5 --fl0_coords0 "22.5,102.1,10500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_SEA_5_${datetime_safe}"
   done

elif [ $loc == "SA" ]; then
   datetimes=("2022-10-01T12:00:00" "2022-01-01T12:00:00")
   t0_fl=("2022-10-01T13:00:00" "2022-01-01T13:00:00")

   for i in "${!datetimes[@]}"; do
      datetime=${datetimes[$i]}
      t0_fl_value=${t0_fl[$i]}

      # Replace colons with underscores
      datetime_safe=$(echo $datetime | tr ':' '_')

      ./run_gpat_sensitivity.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                              --n_ac 1 --fl0_coords0 " -27.5,-67.9,13500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_SA_1_${datetime_safe}"
      ./run_gpat_sensitivity.py --lat_bounds " -28.0,-27.0" --lon_bounds " -68.0,-67.0" --alt_bounds "13000,14000" --t0_sim "$datetime" --t0_fl "$t0_fl_value" \
                              --n_ac 5 --fl0_coords0 " -27.5,-67.9,13500" --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_locs_SA_5_${datetime_safe}"
   done
fi