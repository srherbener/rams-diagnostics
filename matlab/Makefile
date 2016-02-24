#
# make the various script sets
#
#
# Targets
#

all: atex/Makefile atex2/Makefile common/Makefile ideal_etc/Makefile m_map/Makefile metpack/Makefile nctoolbox/Makefile obj_anal/Makefile ts_debby/Makefile nasa_wisc/Makefile tempest/Makefile

%/Makefile: .FORCE
	$(MAKE) -C $(@D)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
