
#Find the nf-config executable from netcdf-fortran that gives us the netcdf setup
NFCONFIG=$(which nf-config)

#Pick a compiler and flags that will work with our netcdf setup, hopefully
FCOMP=$(shell nf-config --fc)
FFLAG=$(shell nf-config --fflags)
FLIBS=$(shell nf-config --flibs)
CLIBS=$(shell nc-config --libs)

#Annoyingly the above does not provide a LD_LIBRARY_PATH type result, so to run things
#one will still need to add them, this should do:

CLPATH=$(shell nc-config --libdir)
FLPATH=$(shell nf-config --prefix)

# We have to guess that nf-config libs are in /lib and not /lib64, nf-config has no --libdir
LPATH=${CLPATH}:${FLPATH}/lib

RUNDIR=/user/home/aesbr/pycontrails_kt/pycontrails

#This is a horrible hack because 'traditional' preprocessing in 30yo fortran code is not fun.
pathsetup.finc: PHONY
	echo "        IERR = NF90_OPEN( '${RUNDIR}/models/files/boxm_ds.nc'  , NF90_WRITE, NCID)" > pathsetup.finc

#Anything relying on PHONY will ALWAYS be rebuilt when make us run
PHONY:

boxm: boxm.f90 Makefile pathsetup.finc
	${FCOMP} -cpp -o boxm boxm.f90 ${FFLAG} ${FLIBS} ${CLIBS} -DPYCONTRAILSDIR=/user/home/aesbr/pycontrails_kt/pycontrails

test: boxm
	LD_LIBRARY_PATH=${LPATH}:${LD_LIBRARY_PATH} ./boxm

#To just see the .f90 file preprocessed run make boxmppre (it'll output to the terminal)
boxmpre: pathsetup.finc
	${FCOMP} -cpp boxm.f90
