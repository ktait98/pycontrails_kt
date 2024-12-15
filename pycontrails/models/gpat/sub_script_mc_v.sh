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

## Params to vary
#   fl_params["n_ac"]: [1, 2, 5, 10]

#   fl_params["sep_dist"][0]: [100, 1000, 2000, 5000, 10000] # dx [m]

#   fl_params["sep_dist"][1]: [0, 100, 1000] # dy [m]

#   plume_params["n_slices"]: [5, 10, 20, 50, 100] # no. of gaussian slices

#   plume_params["max_age"]: [1, 2, 5, 10, 12] # max age of plume waypoints [hours]

#   plume_params["dt_integration"]: [1, 2, 5, 10] # plume simulation int time ["minutes"]

#   plume_params["hres_pl"]: [0.01, 0.02, 0.05, 0.1, 0.5] # plume hres [degrees]
#   plume_params["hres_sim"]: [0.01, 0.02, 0.05, 0.1, 0.5] # chem sim hres [degrees]

# Define the parameter ranges
n_ac=("1" "5" "10")
fl0_heading=("60" "75" "90")
fl0_coords0=("0.225, 0.1, 10500" "0.375, 0.1, 10500" "0.5,0.1,10500")
sep_dist=("0,0,0" "2000,0,0" "5000,0,0" "10000,0,0" "1000,100,0" "1000,1000,0")
n_slices=("20")
max_age=("1" "2" "12")
hres1=("0.025")
hres2=("0.1" "0.5")

# Base case values
base_n_ac="2"
base_fl0_heading="45"
base_fl0_coords0="0.1,0.1,10500"
base_sep_dist="1000,0,0"
base_n_slices="10"
base_max_age="6"
base_hres="0.05"

# Command argument to specify the parameter to vary
param=$1

# Function to run the script with specified parameters
run_script() {
    local n_ac=$1
    local fl0_heading=$2
    local fl0_coords0=$3
    local sep_dist=$4
    local n_slices=$5
    local max_age=$6
    local hres=$7

    IFS=',' read -r -a sep_dist_array <<< "$sep_dist"

   dx=${sep_dist_array[0]}
   dy=${sep_dist_array[1]}
   dz=${sep_dist_array[2]}

    ./run_gpat_mc_v.py --n_ac "$n_ac" --fl0_heading "$fl0_heading" --fl0_coords0 "$fl0_coords0" --sep_dist "$sep_dist" --n_slices "$n_slices" --max_age "$max_age" \
    --hres_sim "$hres" --hres_pl "$hres" --run_gpat --job_id "mc_v_${n_ac}_${fl0_heading}_${dx}_${dy}_${n_slices}_${max_age}_${hres}"
}

# Loop through the specified parameter range
if [ "$param" == "base_case" ]; then
      run_script "$base_n_ac" "$base_fl0_heading" "$base_fl0_coords0" "$base_sep_dist" "$base_n_slices" "$base_max_age" "$base_hres"
elif [ "$param" == "n_ac" ]; then
      for value in "${n_ac[@]}"; do
         run_script "$value" "$base_fl0_heading" "$base_fl0_coords0" "$base_sep_dist" "$base_n_slices" "$base_max_age" "$base_hres"
      done
elif [ "$param" == "fl0_heading" ]; then
    for i in "${!fl0_heading[@]}"; do
        run_script "$base_n_ac" "${fl0_heading[$i]}" "${fl0_coords0[$i]}" "$base_sep_dist" "$base_n_slices" "$base_max_age" "$base_hres"
    done
elif [ "$param" == "sep_dist" ]; then
    for value in "${sep_dist[@]}"; do
        run_script "$base_n_ac" "$base_fl0_heading" "$base_fl0_coords0" "$value" "$base_n_slices" "$base_max_age" "$base_hres"
    done
elif [ "$param" == "n_slices" ]; then
    for value in "${n_slices[@]}"; do
        run_script "$base_n_ac" "$base_fl0_heading" "$base_fl0_coords0" "$base_sep_dist" "$value" "$base_max_age" "$base_hres"
    done
elif [ "$param" == "max_age" ]; then
    for value in "${max_age[@]}"; do
        run_script "$base_n_ac" "$base_fl0_heading" "$base_fl0_coords0" "$base_sep_dist" "$base_n_slices" "$value" "$base_hres"
    done
elif [ "$param" == "hres1" ]; then
    for value in "${hres1[@]}"; do
        run_script "$base_n_ac" "$base_fl0_heading" "$base_fl0_coords0" "$base_sep_dist" "$base_n_slices" "$base_max_age" "$value"
    done
elif [ "$param" == "hres2" ]; then
    for value in "${hres2[@]}"; do
        run_script "$base_n_ac" "$base_fl0_heading" "$base_fl0_coords0" "$base_sep_dist" "$base_n_slices" "$base_max_age" "$value"
    done
else
    echo "Invalid parameter specified. Please choose from: n_ac, fl0_heading, fl0_coords0, sep_dist, n_slices."
fi