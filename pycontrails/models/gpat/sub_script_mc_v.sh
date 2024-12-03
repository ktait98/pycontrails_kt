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
# Steve 20/11/24
#./run_gpat_mc_v.py --n_ac 1 --run_gpat true
# ./run_gpat_mc_v.py --n_ac 5 --run_gpat true
# ./run_gpat_mc_v.py --sep_dist 100,0,0 --run_gpat true
# ./run_gpat_mc_v.py --sep_dist 2000,0,0 --run_gpat true
# ./run_gpat_mc_v.py --sep_dist 1000,1000,0 --run_gpat true
# ./run_gpat_mc_v.py --n_slices 5 --run_gpat true
# ./run_gpat_mc_v.py --n_slices 25 --run_gpat true
# ./run_gpat_mc_v.py --n_slices 50 --run_gpat true
# ./run_gpat_mc_v.py --max_age 12 --run_gpat true
# ./run_gpat_mc_v.py --dt_integration 20  --run_gpat true
# ./run_gpat_mc_v.py --dt_integration 600 --run_gpat true
# ./run_gpat_mc_v.py --hres_pl 0.01 --hres_sim 0.01  --run_gpat true
# ./run_gpat_mc_v.py --hres_pl 0.5  --hres_pl 0.5    --run_gpat true

# Delete mc_v_2_1000_0_10_2_0_0.05 and mc_v_2_1000_0_10_2_10_0.05

# ./run_gpat_mc_v.py --run_gpat true
# ./run_gpat_mc_v.py --dt_integration 1  --run_gpat true
# ./run_gpat_mc_v.py --dt_integration 10 --run_gpat true

./run_gpat_clean.py

