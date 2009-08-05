#
# Config
#

.SUFFIXES: .f90 .o

DEST = $(HOME)/bin
BINS = $(DEST)/azavg \
       $(DEST)/cfad \
       $(DEST)/colint \
       $(DEST)/sfcwind \
       $(DEST)/hslice \
       $(DEST)/vslice

HEADERS = gDataTypes.h

F_COMP = pgf90
#F_COMP_OPTS = -Mupcase -O2 -Ktrap=fp
F_COMP_OPTS = -Mupcase -g -Ktrap=fp
INCLUDES =

LOADER = pgf90
LOADER_OPTS = -v
LIBS =

#
# Targets
#

all: $(BINS)

$(DEST)/azavg: gdata.o azavg_utils.o azavg.o
	$(LOADER) -o $(DEST)/azavg gdata.o azavg_utils.o azavg.o $(LOADER_OPTS) $(LIBS)

$(DEST)/cfad: gdata.o  azavg_utils.o cfad.o
	$(LOADER) -o $(DEST)/cfad gdata.o azavg_utils.o cfad.o $(LOADER_OPTS) $(LIBS)

$(DEST)/colint: gdata.o colint.o
	$(LOADER) -o $(DEST)/colint gdata.o colint.o $(LOADER_OPTS) $(LIBS)

$(DEST)/sfcwind: gdata.o sfcwind.o
	$(LOADER) -o $(DEST)/sfcwind gdata.o sfcwind.o $(LOADER_OPTS) $(LIBS)

$(DEST)/hslice: gdata.o hslice.o
	$(LOADER) -o $(DEST)/hslice gdata.o hslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/vslice: gdata.o vslice.o
	$(LOADER) -o $(DEST)/vslice gdata.o vslice.o $(LOADER_OPTS) $(LIBS)

.f90.o: $(HEADERS)
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

clean:
	rm -f *.o
	rm -f $(BINS)
