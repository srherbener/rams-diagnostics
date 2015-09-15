#
# Config
#

DEST = $(HOME)/etc/matlab/ts_debby
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/TsDebby/SST_OBS

SCRIPTS = \
  $(DEST)/ExtractSampleProfile.m \
  $(DEST)/GenFsAverages.m \
  $(DEST)/GenFsFactors.m \
  $(DEST)/GenMoistDpFiles.m \
  $(DEST)/GenMoistRamsInfo.m \
  $(DEST)/GenMoistRamsInit.m \
  $(DEST)/GenNonSalDpFiles.m \
  $(DEST)/GenNonSalRamsInfo.m \
  $(DEST)/GenNonSalRamsInit.m \
  $(DEST)/GenStormMeasTseries.m \
  $(DEST)/GenStormXsections.m \
  $(DEST)/GenTsdHistMeas.m \
  $(DEST)/LocateSal.m \
  $(DEST)/PlotAeroProfs.m \
  $(DEST)/PlotMaxWindAll.m \
  $(DEST)/PlotMinPressAll.m \
  $(DEST)/PlotTracksAll.m \
  $(DEST)/PlotStormStruct.m \
  $(DEST)/TwoDimAnimation.m \
  $(DEST)/SetPath.m

CONFIGS = \
  $(DIAG_DIR)/FigFsConfig \
  $(DIAG_DIR)/FigPcprateConfig \
  $(DIAG_DIR)/FigVintCondConfig \
  $(DIAG_DIR)/FigHovmollerConfig \
  $(DIAG_DIR)/FigThetaEXsectionConfig \
  $(DIAG_DIR)/FigCoolXsectionConfig \
  $(DIAG_DIR)/FigHeatXsectionConfig \
  $(DIAG_DIR)/FigureConfig

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)

$(DIAG_DIR)/%Config: %Config
	cp $(<) $(@)
