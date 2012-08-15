#
# Config
#

.SUFFIXES: .f90 .mod .o

AZAVG_HOME = $(HOME)/src/diagnostics/tsavg

# Location of results of building the RHDF5 library (ie, running
# make in the 'build' directory of this package
RHDF5_DIR = $(HOME)/etc/fortran
RHDF5_BASE = rhdf5_utils
RHDF5_LIBS = -L$(RHDF5_DIR) -l$(RHDF5_BASE)
RHDF5_INCS = -I$(RHDF5_DIR)

RHDF5_LIB_FILE = $(RHDF5_DIR)/lib$(RHDF5_BASE).a
RHDF5_MOD_FILE = $(RHDF5_DIR)/$(RHDF5_BASE).mod

# diagnostics utilities library
DIAG_DIR = $(HOME)/etc/fortran
DIAG_BASE = diag_utils
DIAG_LIBS = -L$(DIAG_DIR) -l$(DIAG_BASE)
DIAG_INCS = -I$(DIAG_DIR)

DIAG_LIB_FILE = $(DIAG_DIR)/lib$(DIAG_BASE).a
DIAG_MOD_FILE = $(DIAG_DIR)/$(DIAG_BASE).mod

# HDF5 library
HDF5_LIBS=-L$(HOME)/hdf5-1.8.6-macosx64-static/lib -lhdf5 -lhdf5_hl -lsz -lz -lm

# Compiler configuration. These should match what is being set
# in your RAMS and REVU configuration.
INCLUDES = $(RHDF5_INCS) $(DIAG_INCS)

F_COMP = pgf90
F_OPTS = -Munroll -Mnoframe -O2 -pc 64 -tp x64 -Mfree -fastsse -time
F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

LOADER = pgf90
LOADER_OPTS=-v -ldl -lc -lgcc_eh -lpthread -lsz
LOADER_LIBS = $(RHDF5_LIBS) $(HDF5_LIBS) $(DIAG_LIBS)



DEST = $(HOME)/bin
BINS = $(DEST)/tsavg
SRC_DIR = $(AZAVG_HOME)/src

#SCRIPTS = $(DEST)/run_tsavg
SCRIPTS = 

#
# Targets
#

all: $(BINS) $(SCRIPTS)

$(DEST)/tsavg: tsavg.o
	$(LOADER) -o $(DEST)/tsavg tsavg.o $(LOADER_OPTS) $(LOADER_LIBS)

tsavg.o: $(SRC_DIR)/tsavg.f90 $(RHDF5_LIB_FILE) $(RHDF5_MOD_FILE) $(DIAG_LIB_FILE) $(DIAG_MOD_FILE)
	$(F_COMMAND) $(<)

$(DEST)/run_tsavg: $(SRC_DIR)/run_tsavg
	cp $(SRC_DIR)/run_tsavg $(DEST)/run_tsavg

.f90.o:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

.f90.mod:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

#include ./deps.mk

depend:
	makedepf90 -m "%f.mod" -I $(SRC_DIR) -I $(RHDF5_INCS) $(SRC_DIR)/*.f90 > deps.mk

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(BINS)
