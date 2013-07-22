#
# make the various script sets
#
#
# Targets
#

F_UTILS = fortran/utils
F_DIAGS = fortran/azavg fortran/diag_filter fortran/gen_flux fortran/gen_moments fortran/hdata_op fortran/join_hdata fortran/sig_proc
# fortran/tsavg

all: fortran

fortran: $(F_UTILS) $(F_DIAGS)

fortran/%: .FORCE
	$(MAKE) -C $(@)/build

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
