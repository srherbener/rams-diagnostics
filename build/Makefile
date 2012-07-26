#
# Makefile to build rhdf5_utils library and FORTRAN module file
#

# Change RHDF5_HOME depending on where you loaded this package
DIAG_UTILS_HOME = $(HOME)/src/diagnostics/utils

# Install directory should be set to where you want to install
# the resulting library and module files.
DIAG_UTILS_INSTALL_DIR = $(HOME)/etc/fortran

# Location of results of building the RHDF5 library (ie, running
# make in the 'build' directory of this package
RHDF5_DIR = $(HOME)/etc/rhdf5
RHDF5_BASE = rhdf5_utils
RHDF5_LIBS = -L$(RHDF5_DIR) -l$(RHDF5_BASE)
RHDF5_INCS = -I$(RHDF5_DIR)

# Compiler config
INCLUDES = $(RHDF5_INCS)

F_COMP = pgf90
F_OPTS = -Munroll -Mnoframe -O2 -pc 64 -tp x64 -Mfree -fastsse -time
F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

C_COMP = gcc
C_OPTS = -O3 -DUNDERSCORE -DLITTLE -D$(CMACH) $(HDF5_DEFS)
C_COMMAND = $(C_COMP) -c $(C_OPTS) $(INCLUDES)

ARCHIVE = ar rs

# Shouldn't need to change anything below this line.
SRC_DIR = $(DIAG_UTILS_HOME)/src
ARC = libdiag_utils.a
MOD = diag_utils.mod
DEST_ARC = $(DIAG_UTILS_INSTALL_DIR)/$(ARC)
DEST_MOD = $(DIAG_UTILS_INSTALL_DIR)/$(MOD)

OBJS = $(ARC)($(SRC_DIR)/diag_utils.o)

#
# Targets
#

all: $(ARC) $(MOD)

$(ARC): $(OBJS)

install: $(DIAG_UTILS_INSTALL_DIR) $(DEST_ARC) $(DEST_MOD)

$(DIAG_UTILS_INSTALL_DIR):
	mkdir -p $(DIAG_UTILS_INSTALL_DIR)

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

