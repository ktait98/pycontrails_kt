#!/bin/bash

#SBATCH --job-name=BoxModels
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:10:0
#SBATCH --mem=100M
#SBATCH --account=aero004481

#aero004301 
#aero004481 
   #default 
#isys015562


PYCONTRAILSDIR=/user/home/aesbr/pycontrails_kt/pycontrails ./boxm_script.py
