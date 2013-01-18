#
# Config
#

DEST = $(HOME)/etc/matlab/tc_seed
DIAG_DIR = /Volumes/data/RAMS/TcSeeding

SCRIPTS = \
          $(DEST)/AerosolAnimation.m \
          $(DEST)/GenHistOpenSpace.m \
          $(DEST)/GenKeVtPlots.m \
          $(DEST)/MakeDataAvgFig.m \
          $(DEST)/MakeSinglePanelFigs.m \
          $(DEST)/MakeStormEvolFig.m \
          $(DEST)/PlotCfads.m \
          $(DEST)/PlotHmeas.m \
          $(DEST)/PlotHmeas2d.m \
          $(DEST)/PlotHmeas3d.m \
          $(DEST)/PlotHmeasSlices.m \
          $(DEST)/PlotInitialVortex.m \
          $(DEST)/PlotProfs.m \
          $(DEST)/PlotRsProfs.m \
          $(DEST)/PlotSampleCcnProfile.m \
          $(DEST)/PlotSampleStormPhases.m \
          $(DEST)/PlotSampleVt.m \
          $(DEST)/PlotSampleVtSPL.m \
          $(DEST)/PlotSpinUpVtSpl.m \
          $(DEST)/PlotStormPhases.m \
          $(DEST)/PlotStormRegions.m \
          $(DEST)/PlotThetaE.m \
          $(DEST)/PlotTsProfs.m \
          $(DEST)/PlotVt.m \
          $(DEST)/SetPath.m

CONFIGS = \
          $(DIAG_DIR)/SST_26/Init24/DiagConfig \
          $(DIAG_DIR)/SST_26/Init24_SD/DiagConfig \
          $(DIAG_DIR)/SST_26/Init54/DiagConfig

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

$(DIAG_DIR)/SST_26/Init24_SD/DiagConfig: DiagConfig.SST_26.Init24_SD
	cp $(<) $(@)

$(DIAG_DIR)/SST_26/Init54/DiagConfig: DiagConfig.SST_26.Init54
	cp $(<) $(@)
