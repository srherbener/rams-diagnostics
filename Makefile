#
# Config
#

.SUFFIXES: .f90 .o

DEST = $(HOME)/bin
BINS = $(DEST)/azavg

HEADERS = gDataTypes.h
OBJS = gdata.o azavg.o

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

$(DEST)/azavg: $(OBJS)
	$(LOADER) -o $(DEST)/azavg $(OBJS) $(LOADER_OPTS) $(LIBS)

.f90.o: $(HEADERS)
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

clean:
	rm -f $(OBJS)
	rm -f $(BINS)
