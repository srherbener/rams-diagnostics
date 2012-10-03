#
# Config
#

DEST = $(HOME)/etc/matlab/tc_seed

SCRIPTS = $(DEST)/GenKeVtPlots.m \
          $(DEST)/PlotControlVtSpl.m \
          $(DEST)/PlotInitialVortex.m \
          $(DEST)/PlotSampleCcnProfile.m \
          $(DEST)/SetPath.m


#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
