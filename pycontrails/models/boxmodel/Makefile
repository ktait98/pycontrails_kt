
#Find the nf-config executable from netcdf-fortran that gives us the netcdf setup
NFCONFIG=$(which nf-config)

#Pick a compiler and flags that will work with our netcdf setup, hopefully
FCOMP=$(shell nf-config --fc)
FFLAG=$(shell nf-config --fflags)
FLIBS=$(shell nf-config --flibs)


# THESE ARE HACKS FOR BLUECRYSTAL WHICH HAS INCOMPLETE INSTALLATIONS OF NETCDF
#
#  When netcdf-fortran is configured properly, the following command will
#  output the needed info.
#
# nf-config --fflags --flibs
#
# However, currently on bluecrystal we need to fudge the library and header locations
# (much of the problem due to module add not having he correct paths below)
CDFCINST=/software/spack/linux-rocky8-broadwell/gcc-12.3.0/netcdf-c-4.9.2-uujg
CDFFINST=/software/spack/linux-rocky8-broadwell/gcc-12.3.0/netcdf-fortran-4.6.1-fdfq
HACKLIB=-L${CDFCINST}/lib -L${CDFFINST}/lib

LPATH=${CDFCINST}/lib:${CDFFINST}/lib

all: boxm

RUNDIR=/user/home/aesbr/pycontrails_kt/pycontrails

pathsetup.finc: PHONY
	echo "        IERR = NF90_OPEN( '${RUNDIR}/models/files/boxm_ds.nc'  , NF90_WRITE, NCID)" > pathsetup.finc

#Anything relying on PHONY will ALWAYS be rebuilt when make us run
PHONY:

boxm: boxm.f90 Makefile pathsetup.finc
	${FCOMP} -cpp -o boxm boxm.f90 ${FFLAG} ${FLIBS} ${HACKLIB} -DPYCONTRAILSDIR=/user/home/aesbr/pycontrails_kt/pycontrails

test: boxm
	LD_LIBRARY_PATH=${LPATH}:${LD_LIBRARY_PATH} ./boxm

#To just see the .f90 file preprocessed run make boxmppre (it'll output to the terminal)
boxmpre: pathsetup.finc
	${FCOMP} -cpp boxm.f90
