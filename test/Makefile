#
# Makefile to build rhdf5_utils tests
#

# Location of results of building the RHDF5 library (ie, running
# make in the 'build' directory of this package
RHDF5_LIBS = -L/Users/steve/etc/rhdf5 -lrhdf5_utils
RHDF5_INCS = -I/Users/steve/etc/rhdf5

# HDF5 library
HDF5_LIBS=-L/Users/steve/hdf5-1.8.6-macosx64-static/lib -lhdf5 -lsz -lz -lm

# Compiler configuration. These should match what is being set
# in your RAMS and REVU configuration.
INCLUDES = $(RHDF5_INCS)

F_COMP = pgf90
F_OPTS = -Munroll -Mnoframe -O2 -pc 64 -tp x64 -Mfree -fastsse -time
F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

C_COMP = gcc
C_OPTS = -O3 -DUNDERSCORE -DLITTLE $(HDF5_DEFS)
C_COMMAND = $(C_COMP) -c $(C_OPTS) $(INCLUDES)

LOADER = pgf90
LOADER_OPTS=-v -ldl -lc -lgcc_eh -lpthread -lsz
LOADER_LIBS = $(RHDF5_LIBS) $(HDF5_LIBS)

# Executables for testing the RHDF5 library.
#  test_cwrite - Doesn't really test the RHDF5 code, but does give a
#                sense of how fast ~1G of data can be written to disk.
#
#  testh5 - Tests with know array contents so one can see how the data
#           gets organized on disk in the HDF5 file.
#           Also does a timed test for roughly 1G of data.

#
# Targets
#

all: test_cwrite testh5

test_cwrite: test_cwrite.c
	$(C_COMP) -o test_cwrite test_cwrite.c
	rm -f test_cwrite.o

testh5: testh5.f90
	$(F_COMMAND) testh5.f90
	$(LOADER) -o testh5 testh5.o $(LOADER_OPTS) $(LOADER_LIBS)
	rm -f testh5.o

clean:
	rm -f *.o
