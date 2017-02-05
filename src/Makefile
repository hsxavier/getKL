
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
		$(OBJ)/CovMatrix.o $(OBJ)/ParameterList.o
	$(COMP) $(LDHEAL) $(OBJ)/getKL.o $(OBJ)/getKL_aux.o $(OBJ)/fitsfunctions.o $(OBJ)/CovMatrix.o \
			$(OBJ)/Utilities.o $(OBJ)/ParameterList.o -o $@ $(GSLLIBS) $(HEALLIBS)

# Objects:

$(OBJ)/getKL.o: $(SRC)/getKL.cpp $(SRC)/getKL_aux.hpp $(SRC)/Utilities.hpp $(SRC)/fitsfunctions.hpp \
		$(SRC)/CovMatrix.hpp
	$(COMP) -c $(VCONTROL) $(CXXHEAL) $(SRC)/getKL.cpp -o $@ $(CXXOMP)

$(OBJ)/getKL_aux.o: $(SRC)/getKL_aux.cpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/getKL_aux.cpp -o $@

$(OBJ)/CovMatrix.o: $(SRC)/CovMatrix.cpp
	$(COMP) -c $(SRC)/CovMatrix.cpp -o $@

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp
	$(COMP) -c $(SRC)/Utilities.cpp -o $@

$(OBJ)/fitsfunctions.o: $(SRC)/fitsfunctions.cpp $(SRC)/Utilities.hpp
	$(COMP) $(CXXHEAL) $(SETHEALDATA) $(CXXFITS) -c $(SRC)/fitsfunctions.cpp -o $@

$(OBJ)/ParameterList.o: $(SRC)/ParameterList.cpp $(SRC)/ParameterList.hpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/ParameterList.cpp -o $@
