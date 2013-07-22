#
# make the various script sets
#
#
# Targets
#

F_UTILS = utils/build/Makefile
F_DIAGS = azavg/build/Makefile diag_filter/build/Makefile gen_flux/build/Makefile gen_moments/build/Makefile hdata_op/build/Makefile join_hdata/build/Makefile sig_proc/build/Makefile tsavg/build/Makefile

all: fortran

fortran: $(F_UTILS) $(F_DIAGS)

%/Makefile: .FORCE
	$(MAKE) -C $(@D)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
