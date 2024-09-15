#!/bin/bash

#SBATCH --job-name=BoxModels
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:0:10
#SBATCH --mem=100M
#SBATCH --account=aero004481

#aero004301 
#aero004481 
   #default 
#isys015562

./boxm_script.py 
