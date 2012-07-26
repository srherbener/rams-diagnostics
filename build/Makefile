#
# Makefile to build rhdf5_utils library and FORTRAN module file
#

# Change RHDF5_HOME depending on where you loaded this package
RHDF5_HOME = $(HOME)/src/rhdf5

# Install directory should be set to where you want to install
# the resulting library and module files.
RHDF5_INSTALL_DIR = $(HOME)/etc/fortran

# Change HDF5_INCS to the location hdf5 is loaded on your machine.
# This should match what is being set in your RAMS and REVU
# configuration (include.mk.*, etc.).
#
# Uncomment HDF5_DEFS below to use the older (1.6) HDF5 API
#HDF5_DEFS=-DH5_USE_16_API
HDF5_INCS=$(HOME)/hdf5-1.8.6-macosx64-static/include

# Change CMACH to match platform you will compile on. This should
# match the CMACH setting in your RAMS and REVU configuration.
CMACH = PC_LINUX1

# Compiler configuration. These should match what is being set
# in your RAMS and REVU configuration.
INCLUDES = -I$(HDF5_INCS)

F_COMP = pgf90
F_OPTS = -Munroll -Mnoframe -O2 -pc 64 -tp x64 -Mfree -fastsse -time
F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

C_COMP = gcc
C_OPTS = -O3 -DUNDERSCORE -DLITTLE -D$(CMACH) $(HDF5_DEFS)
C_COMMAND = $(C_COMP) -c $(C_OPTS) $(INCLUDES)

ARCHIVE = ar rs

# Shouldn't need to change anything below this line.
SRC_DIR = $(RHDF5_HOME)/src
ARC = librhdf5_utils.a
MOD = rhdf5_utils.mod
DEST_ARC = $(RHDF5_INSTALL_DIR)/$(ARC)
DEST_MOD = $(RHDF5_INSTALL_DIR)/$(MOD)

OBJS = $(ARC)($(SRC_DIR)/rhdf5_f2c.o) \
       $(ARC)($(SRC_DIR)/rhdf5_utils.o)

ARC = librhdf5_utils.a

#
# Targets
#

all: $(ARC) $(MOD)

$(ARC): $(OBJS)

install: $(DEST_ARC) $(DEST_MOD)

$(DEST_ARC): $(ARC)
	cp $(ARC) $(DEST_ARC)

$(DEST_MOD): $(MOD)
	cp $(MOD) $(DEST_MOD)

clean:
	rm -f *.a
	rm -f *.o
	rm -f *.mod

#
# Rules
#

.SUFFIXES: .f90 .c .o .a

# f90 rule.

.f90.a:
	@echo ""
	$(F_COMMAND) $<
	$(ARCHIVE) $@ $(<F:.f90=.o)
	rm -f $(<F:.f90=.o)

# c rule.

.c.a:
	@echo ""
	$(C_COMMAND) $<
	$(ARCHIVE) $@ $(<F:.c=.o)
	rm -f $(<F:.c=.o)

