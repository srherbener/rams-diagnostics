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
  $(DEST)/GenStormProfiles.m \
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
  $(DIAG_DIR)/FigCloudXsectionConfig \
  $(DIAG_DIR)/FigRainXsectionConfig \
  $(DIAG_DIR)/FigDnXsectionConfig \
  $(DIAG_DIR)/FigUpXsectionConfig \
  $(DIAG_DIR)/FigVtXsectionConfig \
  $(DIAG_DIR)/FigVrXsectionConfig \
  $(DIAG_DIR)/FigCloudProfileConfig \
  $(DIAG_DIR)/FigRainProfileConfig \
  $(DIAG_DIR)/FigPrisProfileConfig \
  $(DIAG_DIR)/FigSnowProfileConfig \
  $(DIAG_DIR)/FigAggrProfileConfig \
  $(DIAG_DIR)/FigGraupProfileConfig \
  $(DIAG_DIR)/FigHailProfileConfig \
  $(DIAG_DIR)/FigLhfCoolProfileConfig \
  $(DIAG_DIR)/FigLhfHeatProfileConfig \
  $(DIAG_DIR)/FigLhvCoolProfileConfig \
  $(DIAG_DIR)/FigLhvHeatProfileConfig \
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
