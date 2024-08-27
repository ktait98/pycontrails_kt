# Convert Jupyter to script

Maybe start by converting jupyter notebooks to exectuable scripts ?

https://nbconvert.readthedocs.io/en/latest/usage.html#convert-script

Yeah, could run them as a job as a notebook it seems, but unsure on that wisdom
Perhaps first going with:

nbconvert this.ipynb --out script test_pycontrails.py

To convert the notebook to into a program
alternatively
can import and run in python.

No idea which is best tbh.

# Script to run as start of a cluster job might be something like:

```
#!/bin/sh

git clone https://github.com/ktait98/pycontrails_kt

cd pycontrails_kt

git checkout test

# cd to the fortran directory and rebuild the fortan executables for the platform we're running on here

# FIXME: missing

build_code() {
	cd ./where/the/code/is
	iccorgfortran somecode.f90
}

build_code

# At this point we probably want a way to see what jobs need to be run and how to decide which of them to do next
# But for now would suggest just a testrun verification it works sort of thing

./test_pycontrails_kt.py


# If that all works then we can wrap that script in a slurm submission to test ranges of (hyper)parameters
# Probably best to consider running 10 or 20 in a row on each node we get, depending on run times.

That needs to output to separate directories and name everything appropriately

```
# Scripts/programs to focus on
  - boxm_v_1_ac.py is the latest notebook. It is what currently works with the files it calls, in their current state.
  - flight_gen.py is everything aircraft and plume related. Key output is the emi dataset, which is created once the plumes are aggregated to the grid.
  - boxm.py is the python wrapper around the fortran program boxm.f90. This is its own class and does some of the data wrangling between xarray and netcdf to pass over to fortran. 
  - boxm.f90 is the new implementation of the original box model. It allows for parallelised boxm calcs across a whole region of grid cells. Its outputs match very well with the original boxm, but this validation study needs to be expanded ideally.
  - boxm_orig.for is the old box model code (ew). Caused me great pain over the years staring at and dealing with this code. Oh well, we're basically there with the new version now. Some final validation runs and its good riddance!

# Parameters to vary
```
# met params
  - air temp : 210 -> 280 [k]
  - spec hum : 0.001 -> 0.010 [1]
  - rel hum  : this is set to 0.5. It currently doesn't do anything, but if contrails were modelled later down the line this will become important.
  - e+n wind : 0.0 [m/s] (model currently not capable of modelling eastward and northward wind effects due to no mechanism in BOXM to model mixing between boxes).
  - lag tend : 0.0 [m/s] (lagrangian tendency of air pressure is just a fancy way of saying how much vertical mixing there is between cells. Set to zero for same reason.)

# Actually I don't like this. I think both temp and pressure should be correlated with altitude. I don't see any point in testing 280 K at 40k ft, other than to test the model limits. But the main purpose of this sensitivity study is to explore conditions within the realms of feasibility. Therefore, need to develop simple relation to find avg temp, spec hum and pressure for each altitude level. Almost definitely an answer to this in pycontrails - maybe ISA function or something.

# fl params
  - t0 fl   : e.g. pd.to_datetime("2022-01-20 14:00:00") some sensible time after model initialisation. E.g. if doing a 5 day run, maybe input emission after 24 hours?
  - rt fl   : e.g. pd.Timedelta(minutes=60) set to 60 mins. I think all flights should have a runtime equal to the time it takes for the ac to pass through the box. Maybe this can just be set arbitrarily high and let the lat/lon bounds be what constrains the flight simulation.
  - ts fl   : e.g. pd.Timedelta(minutes=2) determines number of flight waypoints. Note distinction between this and plume waypoints, governed by "dt_int" in plume params.
  - ac type : e.g. A320 I think lets leave this as is for now. It's not a parameter that urgently needs varying, as it covers about half the planes in the sky anyway!
  - fl0 speed : 100 -> 300 [50 m/s]
  - fl0 heading : 0 -> 360 [30 deg] I don't think there is much point in varying heading, as the aircraft is subject to relatively homogeneous conditions. It would only make sense to vary this if it was possible to vary the relative heading between leader and follower flights.
  - fl0 coords : e.g. (47.55, -32.5, 12500) lat, lon, alt [deg, deg, m] again this is relatively arbitrary. We want to capture the greatest physical and chemical response, which means its best to map the entire aircraft track, e.g. N->S, E->W etc. 
  - sep dist : e.g. (5000, 2000, 0) dx, dy, dz [m] this is one of the most crucial params to vary, as this determines the degree of plume overlap between two aircraft. Sorry can't quite think of a good range to observe yet. Just do something arbitrary for now and we can discuss later.
  - n ac    : 0 -> 10 i.e. how many aircraft are in the study. Number of follower aircraft is n_ac - 1.

# plume params
  - dt integration : e.g. pd.Timedelta(minutes=5) i don't know how worthwhile it really is to vary this by much. We need to find a sweetspot where its low enough, but can still handle comp times. 5 mins is good for now.
  - max age : e.g. pd.Timedelta(hours=2) plumes last for 2 - 12 hours supposedly. Maybe we vary between that in 1 hour increments?
  - initial depth and width can stay the same for now.
  - shear   : 0.001 -> 0.1 this is also an important param to play with, as this determines how dispersed the plume becomes over time.

# chem params
  - t0 chem : e.g. pd.to_datetime("2022-01-20 12:00:00") need to have a chat about a sensible range of chem start times.
  - rt chem : e.g. pd.Timedelta(hours=6) I'm thinking to be in line with Anwar and Dick's latest paper, do 5 day runs.
  - ts chem : e.g. pd.Timedelta(seconds=20) don't know if I dare vary this once validation study is done, as doing so resulted in widening discrepancies between boxm1.0 and boxm2.0.
  - lat     : -90 -> +90 [5 deg] e.g. lat_bounds: (-90.0, -85.0)
  - lon     : -180 -> +180 [5 deg] e.g. lon_bounds: (45, 50)
  - level   : 962 (1), 861 (2), 759 (3), 658 (4), 556 (5), 454 (6), 353 (7), 251 (8), 150.5 (9) [hPa] refer to Wasiuk et al. to find pressure altitude ranges for each level. Might as well do the study according to these.
  - hres pl : 0.01 -> 0.1 resolution of plume when aggregated. This is to be higher res than hres chem as the whole point is we need to capture the fine details of the plume expansion, and make sure mass is (somewhat) conserved.
  - hres chem : 0.1 -> 1 [deg] this is highly dependent on mem issues as this seemed to be a problem going too high res here. Maybe not on bc4 tho.
  - vres chem : 500 [deg] can keep this fixed for now.


perhaps (2d params as -thing isn't perfect):

  -region: minlat,maxlat,minlon,maxlon

something like this for two regions to run on?

  -region: 23,24,36,37      33,34,63,64

i.e. run jobs for all possible mixes of them

 ./test_pycontrails -muppets kermit -heights 1m -region 23,24,36,37

and all other combos possible.

That should be easily run across multiple nodes with outputs in different directories.

We can run as many as there is CPU time, if this is my last month or two on bc4 I'm happy to
push it pretty hard on my account :) Seems reasonable, fair use scheduler etc.



