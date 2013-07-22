#
# make the various script sets
#
#
# Targets
#

F_UTILS = fortran/utils
F_DIAGS = azavg.cdir diag_filter.cdir gen_flux.cdir gen_moments.cdir hdata_op.cdir join_hdata.cdir sig_proc.cdir tsavg.cdir

all: fortran

fortran: $(F_UTILS)

fortran/%: .FORCE
	$(MAKE) -C $(@)/build

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
