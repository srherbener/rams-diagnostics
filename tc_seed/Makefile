#
# Config
#

DEST = $(HOME)/etc/matlab/tc_seed

SCRIPTS = $(DEST)/GenKeVtPlots.m \
          $(DEST)/PlotCfads.m \
          $(DEST)/PlotControlVtSpl.m \
          $(DEST)/PlotHmeas.m \
          $(DEST)/PlotHmeas3d.m \
          $(DEST)/PlotHmeasSlices.m \
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
