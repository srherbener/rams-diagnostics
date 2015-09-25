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
  $(DEST)/GenStormMeas.m \
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
  $(DIAG_DIR)/FigMaxWindConfig \
  $(DIAG_DIR)/FigVintCondConfig \
  $(DIAG_DIR)/FigHovmollerConfig \
  $(DIAG_DIR)/FigThetaEXsectionConfig \
  $(DIAG_DIR)/FigCoolXsectionConfig \
  $(DIAG_DIR)/FigHeatXsectionConfig \
  $(DIAG_DIR)/FigCloudXsectionConfig \
  $(DIAG_DIR)/FigRainXsectionConfig \
  $(DIAG_DIR)/FigTcondXsectionConfig \
  $(DIAG_DIR)/FigDnXsectionConfig \
  $(DIAG_DIR)/FigUpXsectionConfig \
  $(DIAG_DIR)/FigVtXsectionConfig \
  $(DIAG_DIR)/FigVrXsectionConfig \
  $(DIAG_DIR)/FigCloudCoreProfileConfig \
  $(DIAG_DIR)/FigRainCoreProfileConfig \
  $(DIAG_DIR)/FigPrisCoreProfileConfig \
  $(DIAG_DIR)/FigSnowCoreProfileConfig \
  $(DIAG_DIR)/FigAggrCoreProfileConfig \
  $(DIAG_DIR)/FigGraupCoreProfileConfig \
  $(DIAG_DIR)/FigHailCoreProfileConfig \
  $(DIAG_DIR)/FigTcondCoreProfileConfig \
  $(DIAG_DIR)/FigLhfCoolCoreProfileConfig \
  $(DIAG_DIR)/FigLhfHeatCoreProfileConfig \
  $(DIAG_DIR)/FigLhvCoolCoreProfileConfig \
  $(DIAG_DIR)/FigLhvHeatCoreProfileConfig \
  $(DIAG_DIR)/FigDnDraftCoreProfileConfig \
  $(DIAG_DIR)/FigUpDraftCoreProfileConfig \
  $(DIAG_DIR)/FigThetaeCoreProfileConfig \
  $(DIAG_DIR)/FigCloudRbandProfileConfig \
  $(DIAG_DIR)/FigRainRbandProfileConfig \
  $(DIAG_DIR)/FigPrisRbandProfileConfig \
  $(DIAG_DIR)/FigSnowRbandProfileConfig \
  $(DIAG_DIR)/FigAggrRbandProfileConfig \
  $(DIAG_DIR)/FigGraupRbandProfileConfig \
  $(DIAG_DIR)/FigHailRbandProfileConfig \
  $(DIAG_DIR)/FigTcondRbandProfileConfig \
  $(DIAG_DIR)/FigLhfCoolRbandProfileConfig \
  $(DIAG_DIR)/FigLhfHeatRbandProfileConfig \
  $(DIAG_DIR)/FigLhvCoolRbandProfileConfig \
  $(DIAG_DIR)/FigLhvHeatRbandProfileConfig \
  $(DIAG_DIR)/FigDnDraftRbandProfileConfig \
  $(DIAG_DIR)/FigUpDraftRbandProfileConfig \
  $(DIAG_DIR)/FigThetaeRbandProfileConfig \
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
