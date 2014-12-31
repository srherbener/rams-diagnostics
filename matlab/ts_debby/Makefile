#
# Config
#

DEST = $(HOME)/etc/matlab/ts_debby
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/TsDebby

SCRIPTS = \
  $(DEST)/LocateSal.m \
  $(DEST)/PlotGrids.m \
  $(DEST)/PlotAeroProfs.m \
  $(DEST)/PlotSampleSoundings.m \
  $(DEST)/GenMoistDpFiles.m \
  $(DEST)/GenMoistRamsInit.m \
  $(DEST)/GenMoistRamsInfo.m \
  $(DEST)/GenStormMeas.m \
  $(DEST)/PlotHydro.m \
  $(DEST)/PlotMaxWind.m \
  $(DEST)/PlotMaxWindAll.m \
  $(DEST)/PlotMinPress.m \
  $(DEST)/PlotMinPressAll.m \
  $(DEST)/PlotPressWind.m \
  $(DEST)/PlotSlp.m \
  $(DEST)/PlotSlpAll.m \
  $(DEST)/PlotSpeedtThetae.m \
  $(DEST)/PlotTracks.m \
  $(DEST)/PlotTracksAll.m \
  $(DEST)/PlotWinds.m \
  $(DEST)/TwoDimAnimation.m \
  $(DEST)/SetPath.m

CONFIGS = \
  $(DIAG_DIR)/SST_CONST/DiagConfig \
  $(DIAG_DIR)/SST_OBS/DiagConfig \
  $(DIAG_DIR)/SST_OBS_UBMIN/DiagConfig

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)

$(DIAG_DIR)/SST_CONST/DiagConfig: DiagConfig.SST_CONST
	cp $(<) $(@)

$(DIAG_DIR)/SST_OBS/DiagConfig: DiagConfig.SST_OBS
	cp $(<) $(@)

$(DIAG_DIR)/SST_OBS_UBMIN/DiagConfig: DiagConfig.SST_OBS_UBMIN
	cp $(<) $(@)
