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

# Set the directory permissions recursively
chmod -R 755 /projects/Impact_of_aviation_on_climate/Kieran2024/

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NCLIBS}:${NFLIBS}
export PYCONTRAILSDIR=/user/work/${USER}/pycontrails_kt/pycontrails/

# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 1 --max_age "ID" --run_gpat --job_id "id_vs_plume_NA_1_1000_0_ID_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 1 --max_age 2 --run_gpat --job_id "id_vs_plume_NA_1_1000_0_2_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 2 --max_age "ID" --run_gpat --job_id "id_vs_plume_NA_2_1000_0_ID_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 2 --max_age 2 --run_gpat --job_id "id_vs_plume_NA_2_1000_0_2_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 3 --max_age "ID" --run_gpat --job_id "id_vs_plume_NA_3_1000_0_ID_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 3 --max_age 2 --run_gpat --job_id "id_vs_plume_NA_3_1000_0_2_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 5 --max_age "ID" --run_gpat --job_id "id_vs_plume_NA_5_1000_0_ID_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 5 --max_age 2 --run_gpat --job_id "id_vs_plume_NA_5_1000_0_2_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 10 --max_age "ID" --run_gpat --job_id "id_vs_plume_NA_10_1000_0_ID_0.01_0.05"
# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" --sep_dist "1000,0,0" --n_ac 10 --max_age 2 --run_gpat --job_id "id_vs_plume_NA_10_1000_0_2_0.01_0.05"

# ./run_gpat_id_vs_plume.py --lat_bounds "47.0,48.0" --lon_bounds " -33.0,-32.0" --alt_bounds "12000,13000" --t0_sim "2022-01-20T12:00:00" \
#                           --n_ac 10 --max_age 2 --sep_dist "2000,0,0" --species_in "NO" --run_gpat --shear "0.05" --job_id "id_vs_plume_NA_10_2000_0_2_0.05_0.05"

# rsync -av --remove-source-files --progress outputs/id_vs_plume_NA_2_1000_0_ID_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
rsync -av --remove-source-files --progress outputs/id_vs_plume_NA_2_1000_0_2_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/

rsync -av --remove-source-files --progress outputs/id_vs_plume_NA_3_1000_0_ID_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
rsync -av --remove-source-files --progress outputs/id_vs_plume_NA_3_1000_0_2_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/

rsync -av --remove-source-files --progress outputs/id_vs_plume_NA_5_1000_0_ID_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
rsync -av --remove-source-files --progress outputs/id_vs_plume_NA_5_1000_0_2_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/

rsync -av --remove-source-files --progress outputs/id_vs_plume_NA_10_1000_0_ID_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/
rsync -av --remove-source-files --progress outputs/id_vs_plume_NA_10_1000_0_2_0.01_0.05 /projects/Impact_of_aviation_on_climate/Kieran2024/outputs/