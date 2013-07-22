#
# make the various script sets
#
#
# Targets
#

F_UTILS = utils/build
F_DIAGS = azavg/build diag_filter/build gen_flux/build gen_moments/build hdata_op/build join_hdata/build sig_proc/build tsavg/build

all: fortran

fortran: $(F_UTILS) $(F_DIAGS)

%/build: .FORCE
	$(MAKE) -C $(@)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
