#
# make the various script sets
#
#
# Targets
#

F_UTILS = utils/build/Makefile getoptions/build/Makefile
F_DIAGS = azavg/build/Makefile create_bins/build/Makefile diag_filter/build/Makefile find_tc_center/build/Makefile gen_flux/build/Makefile gen_moments/build/Makefile gen_xsection/build/Makefile hdata_op/build/Makefile join_hdata/build/Makefile sig_proc/build/Makefile tsavg/build/Makefile gen_adv/build/Makefile gen_vint_terms/build/Makefile

all: fortran

fortran: $(F_UTILS) $(F_DIAGS)

%/Makefile: .FORCE
	$(MAKE) -C $(@D)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
