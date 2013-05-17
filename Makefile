#
# make the various script sets
#
#
# Targets
#

all: atex.mdir common.mdir m_map.mdir metpack.mdir obj_anal.mdir tc_seed.mdir ts_debby.mdir

%.mdir: .FORCE
	$(MAKE) -C $(@)

#
# dummy targets to force action every time
#
.FORCE:

.PHONY: .FORCE
