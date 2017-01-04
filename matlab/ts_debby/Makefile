#
# Config
#

DEST = $(HOME)/etc/matlab/ts_debby
DIAG_DIR = $(HOME)/projects/NsfCollab/Work/SST_OBS

SCRIPTS = \
  $(DEST)/PlotGrids.m \
  $(DEST)/PlotDpCondPrecip.m \
  $(DEST)/PlotDpFig7kmTropTseries.m \
  $(DEST)/PlotDpFigSfcTropTseries.m \
  $(DEST)/PlotDpFigGenRes.m \
  $(DEST)/PlotDpFigHovmollers.m \
  $(DEST)/PlotDpFigStormSeq.m \
  $(DEST)/PlotDpFigDustMassTseries.m \
  $(DEST)/PlotDpFigDustMassTseries_ICCP.m \
  $(DEST)/PlotDpFigDustMassTseries_ICCP_2.m \
  $(DEST)/PlotDpFigDustMassRateTseries.m \
  $(DEST)/PlotDpFigDustRemoval.m \
  $(DEST)/PlotDpFigDustTransport.m \
  $(DEST)/PlotDpFigDustHov.m \
  $(DEST)/PlotDpFigTseries.m \
  $(DEST)/PlotFsFigModelSetup.m \
  $(DEST)/PlotFsFigDustProfiles.m \
  $(DEST)/PlotFsFigCfracProfiles.m \
  $(DEST)/PlotFsFigIke.m \
  $(DEST)/PlotFsFigPcprate.m \
  $(DEST)/PlotFsFigTempPtrack.m \
  $(DEST)/PlotFsFigWindPtrack.m \
  $(DEST)/PlotFsFigVtFactors.m \
  $(DEST)/PlotFsFigCloudRainProfiles.m \
  $(DEST)/PlotFsFigHydroThetaProfiles.m \
  $(DEST)/PlotFsFigAnalysisSetup.m \
  $(DEST)/PlotFsFigIrSimStorm.m \
  $(DEST)/PlotFsFigVisSimStorm.m \
  $(DEST)/PlotFsFigSatSimStorm.m \
  $(DEST)/PlotFsFigObsSimDust.m \
  $(DEST)/PlotFsFigWind.m \
  $(DEST)/PlotFsFigThetaeXsection.m \
  $(DEST)/PlotFsFigThetaXsection.m \
  $(DEST)/PlotFsFigBruntVaisalaXsection.m \
  $(DEST)/PlotFsFigMseXsection.m \
  $(DEST)/PlotFsFigDustPtrackXsection.m \
  $(DEST)/PlotFsFigCondPtrackXsection.m \
  $(DEST)/PlotFsFigDustStrackXsection.m \
  $(DEST)/PlotFsFigCondStrackXsection.m \
  $(DEST)/PlotFsFigVaporXsection.m \
  $(DEST)/PlotFsFigTempWindXsection.m \
  $(DEST)/PlotFsFigTempWindByDustXsection.m \
  $(DEST)/PlotFsFigEvapThetaeXsection.m \
  $(DEST)/PlotFsFigStormStats.m \
  $(DEST)/PlotFsFigStormVt.m \
  $(DEST)/PlotFsFigStormStruct.m \
  $(DEST)/PlotFsFigTseries.m \
  $(DEST)/PlotFsFigBgraph.m \
  $(DEST)/PlotFsFigXsection.m \
  $(DEST)/PlotFsFigHwindVectors.m \
  $(DEST)/PlotFsFigDiagWmap.m \
  $(DEST)/PlotFsFigImageWmap.m \
  $(DEST)/PlotFsFigLine.m \
  $(DEST)/PlotSampleSoundings.m \
  $(DEST)/PlotSampleSoundingsNearStorm.m \
  $(DEST)/ExtractSampleProfile.m \
  $(DEST)/ExtractHorizWinds.m \
  $(DEST)/GenStabilityMeas.m \
  $(DEST)/GenWindShearMeas.m \
  $(DEST)/GenWindShearField.m \
  $(DEST)/GenWindNoVortex.m \
  $(DEST)/GenReduceData.m \
  $(DEST)/GenPtrackMomentum.m \
  $(DEST)/GenPtrackAvgs.m \
  $(DEST)/GenDustAdvection.m \
  $(DEST)/GenTotalMassMeas.m \
  $(DEST)/GenResidualMassMeas.m \
  $(DEST)/GenAzHistSumsCounts.m \
  $(DEST)/GenCintMeas.m \
  $(DEST)/GenTsdVintTseries.m \
  $(DEST)/GenAccumSfcDust.m \
  $(DEST)/GenStormDustFields.m \
  $(DEST)/GenTempAnomalies.m \
  $(DEST)/GenTsavgHists.m \
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
  $(DEST)/GenStormDiffXsections.m \
  $(DEST)/GenStormProfiles.m \
  $(DEST)/GenStormHovmollers.m \
  $(DEST)/GenStrackProfiles.m \
  $(DEST)/GenTsdHistMeas.m \
  $(DEST)/GenTsdHdaMeas.m \
  $(DEST)/GenBuoyancy.m \
  $(DEST)/GenTrackLine.m \
  $(DEST)/ReadGridLoc.m \
  $(DEST)/ReportIntMassNumbers.m \
  $(DEST)/ReportDpMeas.m \
  $(DEST)/ReportDpMidLevelMeas.m \
  $(DEST)/ReportFsMeas.m \
  $(DEST)/RunGenFigures.m \
  $(DEST)/LocateSal.m \
  $(DEST)/PlotAeroProfs.m \
  $(DEST)/PlotDustProfs.m \
  $(DEST)/PlotMaxWindAll.m \
  $(DEST)/PlotMinPressAll.m \
  $(DEST)/PlotTracksAll.m \
  $(DEST)/PlotStormStruct.m \
  $(DEST)/PlotWindVectors.m \
  $(DEST)/TwoDimAnimation.m \
  $(DEST)/SetPath.m

CONFIGS = \
  $(DIAG_DIR)/FigFsConfig \
  $(DIAG_DIR)/FigDpPaperConfig \
  $(DIAG_DIR)/FigStrackCloudConfig \
  $(DIAG_DIR)/FigStrackD1NumConfig \
  $(DIAG_DIR)/FigStrackD2NumConfig \
  $(DIAG_DIR)/FigStrackRainConfig \
  $(DIAG_DIR)/FigStrackRelhumConfig \
  $(DIAG_DIR)/FigStrackSpeedtConfig \
  $(DIAG_DIR)/FigStrackThetaConfig \
  $(DIAG_DIR)/FigStrackThetaeConfig \
  $(DIAG_DIR)/FigStrackVaporConfig \
  $(DIAG_DIR)/FigStrackDustCloudConfig \
  $(DIAG_DIR)/FigStrackDustRainConfig \
  $(DIAG_DIR)/FigStrackDustPrisConfig \
  $(DIAG_DIR)/FigStrackDustSnowConfig \
  $(DIAG_DIR)/FigStrackDustAggrConfig \
  $(DIAG_DIR)/FigStrackDustGraupConfig \
  $(DIAG_DIR)/FigStrackDustHailConfig \
  $(DIAG_DIR)/FigPcprateConfig \
  $(DIAG_DIR)/FigWindShearConfig \
  $(DIAG_DIR)/FigMaxWindConfig \
  $(DIAG_DIR)/FigVintCondConfig \
  $(DIAG_DIR)/FigHovmollerConfig \
  $(DIAG_DIR)/FigThetaeXsectionConfig \
  $(DIAG_DIR)/FigThetaXsectionConfig \
  $(DIAG_DIR)/FigThetaRhoXsectionConfig \
  $(DIAG_DIR)/FigTempcXsectionConfig \
  $(DIAG_DIR)/FigRelhumXsectionConfig \
  $(DIAG_DIR)/FigVaporXsectionConfig \
  $(DIAG_DIR)/FigDustXsectionConfig \
  $(DIAG_DIR)/FigHydroXsectionConfig \
  $(DIAG_DIR)/FigDustHydroXsectionConfig \
  $(DIAG_DIR)/FigMphysXsectionConfig \
  $(DIAG_DIR)/FigCloudNumXsectionConfig \
  $(DIAG_DIR)/FigCloudDiamXsectionConfig \
  $(DIAG_DIR)/FigCloudCondXsectionConfig \
  $(DIAG_DIR)/FigCloudEvapXsectionConfig \
  $(DIAG_DIR)/FigRainNumXsectionConfig \
  $(DIAG_DIR)/FigRainDiamXsectionConfig \
  $(DIAG_DIR)/FigRainCondXsectionConfig \
  $(DIAG_DIR)/FigRainEvapXsectionConfig \
  $(DIAG_DIR)/FigAggrDepXsectionConfig \
  $(DIAG_DIR)/FigAggrSubXsectionConfig \
  $(DIAG_DIR)/FigPrisDepXsectionConfig \
  $(DIAG_DIR)/FigPrisSubXsectionConfig \
  $(DIAG_DIR)/FigSnowDepXsectionConfig \
  $(DIAG_DIR)/FigSnowSubXsectionConfig \
  $(DIAG_DIR)/FigGraupDepXsectionConfig \
  $(DIAG_DIR)/FigGraupSubXsectionConfig \
  $(DIAG_DIR)/FigHailDepXsectionConfig \
  $(DIAG_DIR)/FigHailSubXsectionConfig \
  $(DIAG_DIR)/FigCoolXsectionConfig \
  $(DIAG_DIR)/FigHeatXsectionConfig \
  $(DIAG_DIR)/FigBuoyXsectionConfig \
  $(DIAG_DIR)/FigLiqCondXsectionConfig \
  $(DIAG_DIR)/FigLiqEvapXsectionConfig \
  $(DIAG_DIR)/FigDnXsectionConfig \
  $(DIAG_DIR)/FigUpXsectionConfig \
  $(DIAG_DIR)/FigVtXsectionConfig \
  $(DIAG_DIR)/FigVrXsectionConfig \
  $(DIAG_DIR)/FigHydroProfileConfig \
  $(DIAG_DIR)/FigThetaProfileConfig \
  $(DIAG_DIR)/FigThetaeProfileConfig \
  $(DIAG_DIR)/FigThetaRhoProfileConfig \
  $(DIAG_DIR)/FigUpDraftProfileConfig \
  $(DIAG_DIR)/FigDnDraftCoreProfileConfig \
  $(DIAG_DIR)/FigMphysProfileConfig \
  $(DIAG_DIR)/FigCoolProfileConfig \
  $(DIAG_DIR)/FigHeatProfileConfig \
  $(DIAG_DIR)/FigLiqCondProfileConfig \
  $(DIAG_DIR)/FigLiqEvapProfileConfig \
  $(DIAG_DIR)/FigSpeedtCoreProfileConfig \
  $(DIAG_DIR)/FigTempcCoreProfileConfig \
  $(DIAG_DIR)/FigVaporCoreProfileConfig \
  $(DIAG_DIR)/FigDust1CoreProfileConfig \
  $(DIAG_DIR)/FigDust1NumCoreProfileConfig \
  $(DIAG_DIR)/FigDust2CoreProfileConfig \
  $(DIAG_DIR)/FigDust2NumCoreProfileConfig \
  $(DIAG_DIR)/FigDnDraftRbandProfileConfig \
  $(DIAG_DIR)/FigSpeedtRbandProfileConfig \
  $(DIAG_DIR)/FigTempcRbandProfileConfig \
  $(DIAG_DIR)/FigVaporRbandProfileConfig \
  $(DIAG_DIR)/FigDust1RbandProfileConfig \
  $(DIAG_DIR)/FigDust1NumRbandProfileConfig \
  $(DIAG_DIR)/FigDust2RbandProfileConfig \
  $(DIAG_DIR)/FigDust2NumRbandProfileConfig \
  $(DIAG_DIR)/FigDnDraftEnvProfileConfig \
  $(DIAG_DIR)/FigSpeedtEnvProfileConfig \
  $(DIAG_DIR)/FigTempcEnvProfileConfig \
  $(DIAG_DIR)/FigVaporEnvProfileConfig \
  $(DIAG_DIR)/FigDust1EnvProfileConfig \
  $(DIAG_DIR)/FigDust1NumEnvProfileConfig \
  $(DIAG_DIR)/FigDust2EnvProfileConfig \
  $(DIAG_DIR)/FigDust2NumEnvProfileConfig \
  $(DIAG_DIR)/FigDust1NumSampProfileConfig \
  $(DIAG_DIR)/FigDust2NumSampProfileConfig \
  $(DIAG_DIR)/FigDustProfileConfig \
  $(DIAG_DIR)/FigDustHydroProfileConfig \
  $(DIAG_DIR)/FigCpoolDistConfig \
  $(DIAG_DIR)/FigTseriesConfig \
  $(DIAG_DIR)/FigMidLevelRegenConfig \
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
