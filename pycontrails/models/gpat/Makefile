all: boxm

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

RUNDIR=/user/work/${USER}/pycontrails_kt/pycontrails/

#Anything relying on PHONY will ALWAYS be rebuilt when make us run
PHONY:

# Directories to be created
DIRS = inputs outputs

# Target to create directories
.PHONY: create_dirs
create_dirs:
	mkdir -p $(DIRS)


# Target to download files from Google Drive
.PHONY: download_files
download_files: create_dirs
	
	./download_script.sh 'https://drive.google.com/uc?export=download&id=1-uoRiCyNl_e4o9bbF8ZZuiwYqQcXLfv9' 'inputs/species.nc'
	wget --no-check-certificate 'https://drive.google.com/uc?export=download&id=1gEJQ_XH4RJ6SLJCUs3LLDv8Bn1oiROtk' -O inputs/h2o_concs.nc
	wget --no-check-certificate 'https://drive.google.com/uc?export=download&id=12Qu-8pVXLVb-WH1MphzwJSx9y6xfEtuG' -O inputs/air_temperature.nc

	

# Main target to build the boxm executable
boxm: create_dirs download_files boxm.f90 Makefile
	${FCOMP} -cpp -o boxm boxm.f90 ${FFLAG} ${FLIBS} ${CLIBS} -DPYCONTRAILSDIR=${RUNDIR}

# Clean target to remove the executable and directories
.PHONY: clean
clean:
	rm -f boxm
	rm -rf $(DIRS)

test: boxm
	LD_LIBRARY_PATH=${LPATH}:${LD_LIBRARY_PATH} ./boxm
