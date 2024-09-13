#This should expand to something like /home/username/miniconda3/envs/contrrails/bin/nf-config



#NFCONFIG="${CONDA_PREFIX}/${CONDA_DEFAULT_ENV}/contrails/bin/nf-config"
NFCONFIG=$(which nf-config)

FCOMP=$(shell nf-config --fc)
FFLAG=$(shell nf-config --fflags)
FLIBS=$(shell nf-config --flibs)

#The netcdf-fortran library on bluecrystal phase4 requires -lnetcdf to link
# However, netcdf is not installed with netcdf-fortran and nf-config is not
# setup to return the library location of netcdf.so
# For now, will try hacking it in with this:

CDFCINST=/software/spack/linux-rocky8-broadwell/gcc-12.3.0/netcdf-c-4.9.2-uujg
CDFFINST=/software/spack/linux-rocky8-broadwell/gcc-12.3.0/netcdf-fortran-4.6.1-fdfq
HACKLIB=-L${CDFCINST}/lib -L${CDFFINST}/lib

LPATH=${CDFCINST}/lib:${CDFFINST}/lib

all: boxm

test: boxm
	LD_LIBRARY_PATH=${LPATH}:${LD_LIBRARY_PATH} ./boxm


boxm: boxm.f90 Makefile
	${FCOMP} -o boxm boxm.f90 ${FFLAG} ${FLIBS} ${HACKLIB}
