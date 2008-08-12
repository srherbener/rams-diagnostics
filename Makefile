#
# Config
#

.SUFFIXES: .f90 .o

DEST = $(HOME)/bin
BINS = $(DEST)/azavg

F_COMP = pgf90
F_OPTS =
INCLUDES =

LOADER = pgf90
LOADER_OPTS = -v
LIBS =

#
# Targets
#

all: $(BINS)

$(DEST)/azavg: azavg.o
	$(LOADER) -o $(DEST)/azavg azavg.o $(LOADER_OPTS) $(LIBS)

.f90.o:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

clean:
	rm -f *.o
	rm -f $(BINS)
