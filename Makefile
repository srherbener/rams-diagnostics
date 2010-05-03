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
       $(DEST)/vslice \
       $(DEST)/tr_winds \
       $(DEST)/PDF_pools_allexps

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

$(DEST)/azavg: azavg.o gdata.o azavg_utils.o
	$(LOADER) -o $(DEST)/azavg gdata.o azavg_utils.o azavg.o $(LOADER_OPTS) $(LIBS)

$(DEST)/cfad: cfad.o gdata.o azavg_utils.o
	$(LOADER) -o $(DEST)/cfad gdata.o azavg_utils.o cfad.o $(LOADER_OPTS) $(LIBS)

$(DEST)/colint: colint.o gdata.o
	$(LOADER) -o $(DEST)/colint gdata.o colint.o $(LOADER_OPTS) $(LIBS)

$(DEST)/sfcwind: sfcwind.o gdata.o
	$(LOADER) -o $(DEST)/sfcwind gdata.o sfcwind.o $(LOADER_OPTS) $(LIBS)

$(DEST)/hslice: hslice.o gdata.o
	$(LOADER) -o $(DEST)/hslice gdata.o hslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/vslice: vslice.o gdata.o
	$(LOADER) -o $(DEST)/vslice gdata.o vslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/tr_winds: tr_winds.o gdata.o azavg_utils.o
	$(LOADER) -o $(DEST)/tr_winds gdata.o azavg_utils.o tr_winds.o $(LOADER_OPTS) $(LIBS)

$(DEST)/PDF_pools_allexps: PDF_pools_allexps.o
	$(LOADER) -o $(DEST)/PDF_pools_allexps PDF_pools_allexps.o $(LOADER_OPTS) $(LIBS)

.f90.o: $(HEADERS)
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

clean:
	rm -f *.o
	rm -f $(BINS)
