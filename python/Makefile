#
# make the various script sets
#
#
# Targets
#

all: tutorial/Makefile common/Makefile ts_debby/Makefile

%/Makefile: .FORCE
	$(MAKE) -C $(@D)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
