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



n_ac=("2" "3" "10")
sep_dists=("0" "100" "500" "1000" "2000" "3000" "4000" "5000" "10000")

sep_dist=$1

if [ $sep_dist == "5_ac_10000" ]; then
   ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac 5 --sep_dist "10000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_5_10000_0_6_0.01_0.05"
fi

if [ $sep_dist == "0" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "0,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_0_0_6_0.01_0.05"
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_0_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_0_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
elif [ $sep_dist == "100" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "100,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_100_0_6_0.01_0.05"
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_100_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_100_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
elif [ $sep_dist == "500" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "500,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_500_0_6_0.01_0.05"
   #    rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_500_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
   #    rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_500_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
elif [ $sep_dist == "1000" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "1000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_1000_0_6_0.01_0.05"
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_1000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_1000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
elif [ $sep_dist == "2000" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "2000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_2000_0_6_0.01_0.05"
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_2000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_2000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
elif [ $sep_dist == "3000" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "3000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_3000_0_6_0.01_0.05"
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_3000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_3000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
elif [ $sep_dist == "4000" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "4000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_4000_0_6_0.01_0.05"
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_4000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_4000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
elif [ $sep_dist == "5000" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "5000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_5000_0_6_0.01_0.05"
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_5000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_5000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
elif [ $sep_dist == "10000" ]; then
   for n_ac in "${n_ac[@]}"; do
      ./run_gpat_sensitivity.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
                           --n_ac $n_ac --sep_dist "10000,0,0" --max_age 6 --run_gpat --job_id "sensitivity_NA_${n_ac}_10000_0_6_0.01_0.05"
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/inputs/sensitivity_NA_${n_ac}_10000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/inputs/
      # rsync -av --remove-source-files --progress /user/work/${USER}/pycontrails_kt/pycontrails/models/gpat/outputs/sensitivity_NA_${n_ac}_10000_0_6_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
   done
fi