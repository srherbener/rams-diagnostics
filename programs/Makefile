#
# Config
#

.SUFFIXES: .f90 .mod .o

DEST = $(HOME)/bin
BINS = $(DEST)/azavg \
       $(DEST)/cfad \
       $(DEST)/colint \
       $(DEST)/densf \
       $(DEST)/ghist \
       $(DEST)/gfilter \
       $(DEST)/sfcwind \
       $(DEST)/hslice \
       $(DEST)/vslice \
       $(DEST)/tr_winds \
       $(DEST)/PDF_pools_allexps \
       $(DEST)/prueba39 \
       $(DEST)/join_gdata \
       $(DEST)/double_pts \
       $(DEST)/mkTemplate \
       $(DEST)/tseries_avg

F_COMP = pgf90
#F_COMP_OPTS = -tp px-64 -pc 64 -Mupcase -O2 -Ktrap=fp
F_COMP_OPTS = -tp px-64 -pc 64 -Mupcase -g -Ktrap=fp
INCLUDES =

LOADER = pgf90
LOADER_OPTS = -v
LIBS =

#
# Targets
#

all: $(BINS)

$(DEST)/azavg: azavg.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/azavg gdata_utils.o azavg_utils.o azavg.o $(LOADER_OPTS) $(LIBS)

$(DEST)/cfad: cfad.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/cfad gdata_utils.o azavg_utils.o cfad.o $(LOADER_OPTS) $(LIBS)

$(DEST)/colint: colint.o gdata_utils.o
	$(LOADER) -o $(DEST)/colint gdata_utils.o colint.o $(LOADER_OPTS) $(LIBS)

$(DEST)/densf: densf.o gdata_utils.o
	$(LOADER) -o $(DEST)/densf gdata_utils.o densf.o $(LOADER_OPTS) $(LIBS)

$(DEST)/ghist: ghist.o gdata_utils.o
	$(LOADER) -o $(DEST)/ghist gdata_utils.o ghist.o $(LOADER_OPTS) $(LIBS)

$(DEST)/gfilter: gfilter.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/gfilter gdata_utils.o azavg_utils.o gfilter.o $(LOADER_OPTS) $(LIBS)

$(DEST)/sfcwind: sfcwind.o gdata_utils.o
	$(LOADER) -o $(DEST)/sfcwind gdata_utils.o sfcwind.o $(LOADER_OPTS) $(LIBS)

$(DEST)/hslice: hslice.o gdata_utils.o
	$(LOADER) -o $(DEST)/hslice gdata_utils.o hslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/vslice: vslice.o gdata_utils.o
	$(LOADER) -o $(DEST)/vslice gdata_utils.o vslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/tr_winds: tr_winds.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/tr_winds gdata_utils.o azavg_utils.o tr_winds.o $(LOADER_OPTS) $(LIBS)

$(DEST)/PDF_pools_allexps: PDF_pools_allexps.o gdata_utils.o
	$(LOADER) -o $(DEST)/PDF_pools_allexps gdata_utils.o PDF_pools_allexps.o $(LOADER_OPTS) $(LIBS)

$(DEST)/prueba39: prueba39.o
	$(LOADER) -o $(DEST)/prueba39 prueba39.o $(LOADER_OPTS) $(LIBS)

$(DEST)/join_gdata: join_gdata.o gdata_utils.o
	$(LOADER) -o $(DEST)/join_gdata gdata_utils.o join_gdata.o $(LOADER_OPTS) $(LIBS)

$(DEST)/double_pts: double_pts.o gdata_utils.o
	$(LOADER) -o $(DEST)/double_pts gdata_utils.o double_pts.o $(LOADER_OPTS) $(LIBS)

$(DEST)/tseries_avg: tseries_avg.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/tseries_avg gdata_utils.o azavg_utils.o tseries_avg.o $(LOADER_OPTS) $(LIBS)

$(DEST)/mkTemplate: mkTemplate.o gdata_utils.o
	$(LOADER) -o $(DEST)/mkTemplate gdata_utils.o mkTemplate.o $(LOADER_OPTS) $(LIBS)

.f90.o:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

.f90.mod:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

azavg.o: azavg.f90 gdata_utils.mod azavg_utils.mod
cfad.o: cfad.f90 gdata_utils.mod azavg_utils.mod
colint.o: colint.f90 gdata_utils.mod
densf.o: densf.f90 gdata_utils.mod
double_pts.o: double_pts.f90 gdata_utils.mod
gfilter.o: gfilter.f90 gdata_utils.mod azavg_utils.mod
ghist.o: ghist.f90 gdata_utils.mod
hslice.o: hslice.f90 gdata_utils.mod
join_gdata.o: join_gdata.f90 gdata_utils.mod
mkTemplate.o: mkTemplate.f90 gdata_utils.mod
PDF_pools_allexps.o: PDF_pools_allexps.f90 gdata_utils.mod
sfcwind.o: sfcwind.f90 gdata_utils.mod
tr_winds.o: tr_winds.f90 gdata_utils.mod azavg_utils.mod
tseries_avg.o: tseries_avg.f90 gdata_utils.mod azavg_utils.mod
vslice.o: vslice.f90 gdata_utils.mod

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(BINS)
