#
# make the various script sets
#
#
# Targets
#

all: fortran

fortran: .FORCE
	$(MAKE) -C $(@)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
