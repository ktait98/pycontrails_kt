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
#

CDFCINST=/software/spack/linux-rocky8-broadwell/gcc-12.3.0/netcdf-c-4.9.2-uujg
CDFFINST=/software/spack/linux-rocky8-broadwell/gcc-12.3.0/netcdf-fortran-4.6.1-fdfq
HACKLIB=-L${CDFCINST}/lib -L${CDFFINST}/lib

LPATH=${CDFCINST}/lib:${CDFFINST}/lib

# HACKLIB does seem to work, but it's an awful hack
# and ... it leaves boxm a bit like this:
#
# [aesbr@bc4login3 boxmodel]$ ldd boxm
        # linux-vdso.so.1 (0x00007fe1ad46b000)
        # libnetcdff.so.7 => not found
        # libnetcdf.so.19 => not found
        # libdl.so.2 => /lib64/libdl.so.2 (0x00007fe1ad03b000)
        # libgfortran.so.5 => /software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-vpim/lib64/libgfortran.so.5 (0x00007fe1acd67000)
        # libm.so.6 => /lib64/libm.so.6 (0x00007fe1ac9e5000)
        # libgcc_s.so.1 => /software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-vpim/lib64/libgcc_s.so.1 (0x00007fe1ad430000)
        # libquadmath.so.0 => /software/spack/linux-rocky8-broadwell/gcc-12.3.0/gcc-12.3.0-vpim/lib64/libquadmath.so.0 (0x00007fe1ad3e9000)
        # libc.so.6 => /lib64/libc.so.6 (0x00007fe1ac620000)
        # /lib64/ld-linux-x86-64.so.2 (0x00007fe1ad23f000)
# [aesbr@bc4login3 boxmodel]$ 

all: boxm

test: boxm
	LD_LIBRARY_PATH=${LPATH}:${LD_LIBRARY_PATH} ./boxm


boxm: boxm.f90 Makefile
	#gfortran -o boxm boxm.f90 $(/home/ktait98/miniconda3/envs/contrails/bin/nf-config --fflags) $(/home/ktait98/miniconda3/envs/contrails/bin/nf-config --flibs)
	${FCOMP} -o boxm boxm.f90 ${FFLAG} ${FLIBS} ${HACKLIB}
