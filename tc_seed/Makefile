#
# Config
#

DEST = $(HOME)/etc/matlab/tc_seed
DIAG_DIR = /Volumes/data/RAMS/TcSeeding

SCRIPTS = $(DEST)/GenKeVtPlots.m \
          $(DEST)/PlotCfads.m \
          $(DEST)/PlotControlVtSpl.m \
          $(DEST)/PlotHmeas.m \
          $(DEST)/PlotHmeas2d.m \
          $(DEST)/PlotHmeas3d.m \
          $(DEST)/PlotHmeasSlices.m \
          $(DEST)/PlotInitialVortex.m \
          $(DEST)/PlotProfs.m \
          $(DEST)/PlotSampleCcnProfile.m \
          $(DEST)/PlotThetaE.m \
          $(DEST)/PlotTsProfs.m \
          $(DEST)/PlotVt.m \
          $(DEST)/SetPath.m

CONFIGS = $(DIAG_DIR)/SST_26/Init24/DiagConfig

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)

$(DIAG_DIR)/SST_26/Init24/DiagConfig: DiagConfig.SST_26.Init24
	cp $(<) $(@)
