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

$(DEST)/azavg: azavg.o gdata_mod.o azavg_utils.o
	$(LOADER) -o $(DEST)/azavg gdata_mod.o azavg_utils.o azavg.o $(LOADER_OPTS) $(LIBS)

$(DEST)/cfad: cfad.o gdata_mod.o azavg_utils.o
	$(LOADER) -o $(DEST)/cfad gdata_mod.o azavg_utils.o cfad.o $(LOADER_OPTS) $(LIBS)

$(DEST)/colint: colint.o gdata_mod.o
	$(LOADER) -o $(DEST)/colint gdata_mod.o colint.o $(LOADER_OPTS) $(LIBS)

$(DEST)/densf: densf.o gdata_mod.o
	$(LOADER) -o $(DEST)/densf gdata_mod.o densf.o $(LOADER_OPTS) $(LIBS)

$(DEST)/ghist: ghist.o gdata_mod.o
	$(LOADER) -o $(DEST)/ghist gdata_mod.o ghist.o $(LOADER_OPTS) $(LIBS)

$(DEST)/gfilter: gfilter.o gdata_mod.o azavg_utils.o
	$(LOADER) -o $(DEST)/gfilter gdata_mod.o azavg_utils.o gfilter.o $(LOADER_OPTS) $(LIBS)

$(DEST)/sfcwind: sfcwind.o gdata_mod.o
	$(LOADER) -o $(DEST)/sfcwind gdata_mod.o sfcwind.o $(LOADER_OPTS) $(LIBS)

$(DEST)/hslice: hslice.o gdata_mod.o
	$(LOADER) -o $(DEST)/hslice gdata_mod.o hslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/vslice: vslice.o gdata_mod.o
	$(LOADER) -o $(DEST)/vslice gdata_mod.o vslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/tr_winds: tr_winds.o gdata_mod.o azavg_utils.o
	$(LOADER) -o $(DEST)/tr_winds gdata_mod.o azavg_utils.o tr_winds.o $(LOADER_OPTS) $(LIBS)

$(DEST)/PDF_pools_allexps: PDF_pools_allexps.o gdata_mod.o
	$(LOADER) -o $(DEST)/PDF_pools_allexps gdata_mod.o PDF_pools_allexps.o $(LOADER_OPTS) $(LIBS)

$(DEST)/prueba39: prueba39.o
	$(LOADER) -o $(DEST)/prueba39 prueba39.o $(LOADER_OPTS) $(LIBS)

$(DEST)/join_gdata: join_gdata.o gdata_mod.o
	$(LOADER) -o $(DEST)/join_gdata gdata_mod.o join_gdata.o $(LOADER_OPTS) $(LIBS)

$(DEST)/double_pts: double_pts.o gdata_mod.o
	$(LOADER) -o $(DEST)/double_pts gdata_mod.o double_pts.o $(LOADER_OPTS) $(LIBS)

$(DEST)/tseries_avg: tseries_avg.o gdata_mod.o azavg_utils.o
	$(LOADER) -o $(DEST)/tseries_avg gdata_mod.o azavg_utils.o tseries_avg.o $(LOADER_OPTS) $(LIBS)

$(DEST)/mkTemplate: mkTemplate.o gdata_mod.o
	$(LOADER) -o $(DEST)/mkTemplate gdata_mod.o mkTemplate.o $(LOADER_OPTS) $(LIBS)

.f90.o:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

.f90.mod:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

azavg.o: azavg.f90 gdata_mod.mod
azavg_utils.o: azavg_utils.f90 gdata_mod.mod
cfad.o: cfad.f90 gdata_mod.mod
colint.o: colint.f90 gdata_mod.mod
densf.o: densf.f90 gdata_mod.mod
double_pts.o: double_pts.f90 gdata_mod.mod
gfilter.o: gfilter.f90 gdata_mod.mod
ghist.o: ghist.f90 gdata_mod.mod
hslice.o: hslice.f90 gdata_mod.mod
join_gdata.o: join_gdata.f90 gdata_mod.mod
mkTemplate.o: mkTemplate.f90 gdata_mod.mod
PDF_pools_allexps.o: PDF_pools_allexps.f90 gdata_mod.mod
sfcwind.o: sfcwind.f90 gdata_mod.mod
tr_winds.o: tr_winds.f90 gdata_mod.mod
tseries_avg.o: tseries_avg.f90 gdata_mod.mod
vslice.o: vslice.f90 gdata_mod.mod

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(BINS)
