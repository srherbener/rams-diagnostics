#
# make the various script sets
#
#
# Targets
#

F_UTILS = fortran/utils
F_DIAGS = fortran/azavg fortran/diag_filter
# gen_flux.cdir gen_moments.cdir hdata_op.cdir join_hdata.cdir sig_proc.cdir tsavg.cdir

all: fortran

fortran: $(F_UTILS) $(F_DIAGS)

fortran/%: .FORCE
	$(MAKE) -C $(@)/build

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
