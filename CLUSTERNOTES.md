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

# Paramaters to vary ?

# Making these up for now, example placeholders, but this might be a nice easy text file format way to do it:

space separated variables options, so I guess no spaces in the variable might make this easier.

  -muppets : kermit piggy beaker
  -heights : 1m 20m 20km
  -lat     : 23-27 33-90 4-5
  -lon     : ...

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



