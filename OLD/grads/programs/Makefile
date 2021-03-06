#
# Config
#

.SUFFIXES: .f90 .mod .o

DEST = $(HOME)/bin
BINS = $(DEST)/azavg_grads \
       $(DEST)/cfad_grads \
       $(DEST)/colint_grads \
       $(DEST)/densf \
       $(DEST)/ghist \
       $(DEST)/gfilter \
       $(DEST)/sfcwind_grads \
       $(DEST)/hslice_grads \
       $(DEST)/vslice \
       $(DEST)/tr_winds_grads \
       $(DEST)/PDF_pools_allexps \
       $(DEST)/prueba39 \
       $(DEST)/join_gdata \
       $(DEST)/double_pts \
       $(DEST)/mkTemplate \
       $(DEST)/tseries_avg_grads

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

$(DEST)/azavg_grads: azavg.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/azavg_grads gdata_utils.o azavg_utils.o azavg.o $(LOADER_OPTS) $(LIBS)

$(DEST)/cfad_grads: cfad.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/cfad_grads gdata_utils.o azavg_utils.o cfad.o $(LOADER_OPTS) $(LIBS)

$(DEST)/colint_grads: colint.o gdata_utils.o
	$(LOADER) -o $(DEST)/colint_grads gdata_utils.o colint.o $(LOADER_OPTS) $(LIBS)

$(DEST)/densf: densf.o gdata_utils.o
	$(LOADER) -o $(DEST)/densf gdata_utils.o densf.o $(LOADER_OPTS) $(LIBS)

$(DEST)/ghist: ghist.o gdata_utils.o
	$(LOADER) -o $(DEST)/ghist gdata_utils.o ghist.o $(LOADER_OPTS) $(LIBS)

$(DEST)/gfilter: gfilter.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/gfilter gdata_utils.o azavg_utils.o gfilter.o $(LOADER_OPTS) $(LIBS)

$(DEST)/sfcwind_grads: sfcwind.o gdata_utils.o
	$(LOADER) -o $(DEST)/sfcwind_grads gdata_utils.o sfcwind.o $(LOADER_OPTS) $(LIBS)

$(DEST)/hslice_grads: hslice.o gdata_utils.o
	$(LOADER) -o $(DEST)/hslice_grads gdata_utils.o hslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/vslice: vslice.o gdata_utils.o
	$(LOADER) -o $(DEST)/vslice gdata_utils.o vslice.o $(LOADER_OPTS) $(LIBS)

$(DEST)/tr_winds_grads: tr_winds.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/tr_winds_grads gdata_utils.o azavg_utils.o tr_winds.o $(LOADER_OPTS) $(LIBS)

$(DEST)/PDF_pools_allexps: PDF_pools_allexps.o gdata_utils.o
	$(LOADER) -o $(DEST)/PDF_pools_allexps gdata_utils.o PDF_pools_allexps.o $(LOADER_OPTS) $(LIBS)

$(DEST)/prueba39: prueba39.o
	$(LOADER) -o $(DEST)/prueba39 prueba39.o $(LOADER_OPTS) $(LIBS)

$(DEST)/join_gdata: join_gdata.o gdata_utils.o
	$(LOADER) -o $(DEST)/join_gdata gdata_utils.o join_gdata.o $(LOADER_OPTS) $(LIBS)

$(DEST)/double_pts: double_pts.o gdata_utils.o
	$(LOADER) -o $(DEST)/double_pts gdata_utils.o double_pts.o $(LOADER_OPTS) $(LIBS)

$(DEST)/tseries_avg_grads: tseries_avg.o gdata_utils.o azavg_utils.o
	$(LOADER) -o $(DEST)/tseries_avg_grads gdata_utils.o azavg_utils.o tseries_avg.o $(LOADER_OPTS) $(LIBS)

$(DEST)/mkTemplate: mkTemplate.o gdata_utils.o
	$(LOADER) -o $(DEST)/mkTemplate gdata_utils.o mkTemplate.o $(LOADER_OPTS) $(LIBS)

.f90.o:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

.f90.mod:
	$(F_COMP) $(F_COMP_OPTS) $(INCLUDES) -c $(<)

include ./deps.mk

depend:
	makedepf90 -m "%f.mod" *.f90 > deps.mk

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(BINS)
