
############################
### Variable definitions ###
############################

# Option for debugging with valgrind:
DEBUG = -g -O0
# COMP = g++ $(DEBUG)

# Compiler and flask directory structure:
COMP = g++
BIN  = ../bin
OBJ  = ../obj
SRC  = ../src

# Healpix location:
HEALDIR  =  /home/skems/prog/Healpix_3.31
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib

# Other compiler flags:
CXXOMP   = -fopenmp
GSLLIBS  = -lgsl -lgslcblas
HEALLIBS = -lhealpix_cxx -lcxxsupport -lcfitsio -lsharp -lfftpack -lc_utils -lgomp
FFTWLIBS = -lfftw3 -lm

# Options for Splinter computer:
ifeq ($(HOSTNAME), splinter-login.local)
HEALDIR  = /home/hsxavier/Healpix_3.11
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib
CXXFITS  = -I/usr/include/cfitsio
endif

# Options for Alphacrucis computer:
ifeq ($(HOSTNAME), alphacrucis)
HEALDIR  = /sto/home/hsxavier/Healpix_3.20
HEALDATA = $(HEALDIR)/data
CXXHEAL  = -I$(HEALDIR)/src/cxx/generic_gcc/include
LDHEAL   = -L$(HEALDIR)/src/cxx/generic_gcc/lib
CXXFITS  = -I/sto/home/hsxavier/cfitsio
LDFITS   = -L/sto/home/hsxavier/cfitsio
endif

# Get code version:
COMMIT := $(shell git rev-parse --short HEAD)
ifeq ($(COMMIT), )
COMMIT = missing
endif

# Dynamic defines (#define that are set during compilation time, by passing instruction to compiler):
VCONTROL    = -D COMMIT=\"$(COMMIT)\"
SETHEALDATA = -D HEALPIX_DATA=\"$(HEALDATA)\"



#############################
### Compiling intructions ###
#############################

# General instructions:

all: $(BIN)/getKL

clean:
	rm -f $(BIN)/*
	rm -f $(OBJ)/*


# Executables:

$(BIN)/getKL: $(OBJ)/getKL.o $(OBJ)/getKL_aux.o $(OBJ)/fitsfunctions.o $(OBJ)/Utilities.o \
		$(OBJ)/CovMatrix.o $(OBJ)/ParameterList.o $(OBJ)/interpol.o $(OBJ)/HealpixMapProcess.o \
		$(OBJ)/Cosmology.o $(OBJ)/Integral.o $(OBJ)/FFTrelated.o
	$(COMP) $(LDHEAL) $(OBJ)/getKL.o $(OBJ)/getKL_aux.o $(OBJ)/fitsfunctions.o $(OBJ)/CovMatrix.o \
			$(OBJ)/Utilities.o $(OBJ)/ParameterList.o $(OBJ)/interpol.o $(OBJ)/HealpixMapProcess.o \
		$(OBJ)/Cosmology.o $(OBJ)/Integral.o $(OBJ)/FFTrelated.o -o $@ $(GSLLIBS) $(HEALLIBS) $(FFTWLIBS)


# Objects:

$(OBJ)/getKL.o: $(SRC)/getKL.cpp $(SRC)/getKL_aux.hpp $(SRC)/Utilities.hpp $(SRC)/fitsfunctions.hpp \
		$(SRC)/CovMatrix.hpp $(SRC)/definitions.h $(SRC)/HealpixMapProcess.hpp $(SRC)/Cosmology.hpp \
		$(SRC)/FFTrelated.hpp $(SRC)/output.hpp
	$(COMP) -c $(VCONTROL) $(CXXHEAL) $(SRC)/getKL.cpp -o $@ $(CXXOMP)

$(OBJ)/getKL_aux.o: $(SRC)/getKL_aux.cpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/getKL_aux.cpp -o $@

$(OBJ)/FFTrelated.o: $(SRC)/FFTrelated.cpp $(SRC)/Utilities.hpp $(SRC)/interpol.h
	$(COMP) -c $(SRC)/FFTrelated.cpp -o $@

$(OBJ)/CovMatrix.o: $(SRC)/CovMatrix.cpp
	$(COMP) -c $(SRC)/CovMatrix.cpp -o $@

$(OBJ)/HealpixMapProcess.o: $(SRC)/HealpixMapProcess.cpp $(SRC)/definitions.h
	$(COMP) $(CXXHEAL) -c $(SRC)/HealpixMapProcess.cpp -o $@

$(OBJ)/Cosmology.o: $(SRC)/Cosmology.cpp $(SRC)/Cosmology.hpp $(SRC)/Utilities.hpp $(SRC)/Integral.hpp \
	$(SRC)/ParameterList.hpp $(SRC)/interpol.h
	$(COMP) -c $(SRC)/Cosmology.cpp -o $@

$(OBJ)/Integral.o: $(SRC)/Integral.cpp $(SRC)/Integral.hpp $(SRC)/Utilities.hpp $(SRC)/Cosmology.hpp
	$(COMP) -c $(SRC)/Integral.cpp -o $@

$(OBJ)/interpol.o: $(SRC)/interpol.c
	$(COMP) -c $(SRC)/interpol.c -o $@

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp
	$(COMP) -c $(SRC)/Utilities.cpp -o $@

$(OBJ)/fitsfunctions.o: $(SRC)/fitsfunctions.cpp $(SRC)/Utilities.hpp
	$(COMP) $(CXXHEAL) $(SETHEALDATA) $(CXXFITS) -c $(SRC)/fitsfunctions.cpp -o $@

$(OBJ)/ParameterList.o: $(SRC)/ParameterList.cpp $(SRC)/ParameterList.hpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/ParameterList.cpp -o $@
