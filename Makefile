#
# make the various script sets
#
#
# Targets
#

UTILS_DIR = utils.cdir
CODE_DIRS = azavg.cdir diag_filter.cdir gen_flux.cdir gen_moments.cdir hdata_op.cdir join_hdata.cdir sig_proc.cdir tsavg.cdir

all: $(UTILS_DIR) $(CODE_DIRS)

%.cdir: .FORCE
	$(MAKE) -C $(@)/build

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
