#!/bin/bash

#SBATCH --job-name=BoxModels
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:0:0
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

NCLIBS=`nc-config --libdir`
NFLIBS=`nf-config --prefix`/lib

LD_LIBRARY_PATH=${NCLIBS}:${NFLIBS} PYCONTRAILSDIR=/user/work/${USER}/pycontrails_kt/pycontrails/ ./run_gpat.py
