#
# make the various script sets
#
#
# Targets
#

all: fortran/Makefile grads/Makefile idl/Makefile

%/Makefile: .FORCE
	$(MAKE) -C $(@D)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
