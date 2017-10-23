#
# make the various script sets
#
#
# Targets
#

all: rhdf5_install fortran/Makefile grads/Makefile idl/Makefile matlab/Makefile perl/Makefile python/Makefile utils/Makefile

rhdf5_install: .FORCE
	$(MAKE) -C rhdf5/build install

%/Makefile: .FORCE
	$(MAKE) -C $(@D)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
