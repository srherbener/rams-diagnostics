#
# Config
#

DEST = $(HOME)/etc/matlab/atex
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/MarinePrecip/ATEX

SCRIPTS = $(DEST)/ExtractSample.m \
          $(DEST)/CombineStCloudTypes.m \
          $(DEST)/ReportAvgCtypeCfrac.m \
          $(DEST)/ReportEnsLtss.m \
          $(DEST)/GenEnsembleData.m \
          $(DEST)/GenCloudFracPlots.m \
          $(DEST)/GenFigures.m \
          $(DEST)/GenLegendPlots.m \
          $(DEST)/GenAvgCloudDepth.m \
          $(DEST)/GenAvgLwpToCdepthFilesCtype.m \
          $(DEST)/GenAvgFilesCtype.m \
          $(DEST)/GenAvgLCL.m \
          $(DEST)/GenAlbedo.m \
          $(DEST)/GenBLcouplingData.m \
          $(DEST)/GenBuoyancyFlux.m \
          $(DEST)/GenBlStats.m \
          $(DEST)/GenCloudStructPlots.m \
          $(DEST)/GenDistPlots.m \
          $(DEST)/GenDropSizeDists.m \
          $(DEST)/GenDomAvgCtypeTseries.m \
          $(DEST)/GenHfracCloud.m \
          $(DEST)/GenHdaTseriesFiles.m \
          $(DEST)/GenBarGraphData.m \
          $(DEST)/GenBarGraphFilesCtype.m \
          $(DEST)/GenBarGraphsCtype.m \
          $(DEST)/GenMomentData.m \
          $(DEST)/GenProfileData.m \
          $(DEST)/GenHistData.m \
          $(DEST)/GenHistFiles.m \
          $(DEST)/GenHistFilesCtype.m \
          $(DEST)/GenPdfFiles.m \
          $(DEST)/GenPdfFilesCtype.m \
          $(DEST)/GenPdfPlots.m \
          $(DEST)/GenTurbData.m \
          $(DEST)/GenWcfads.m \
          $(DEST)/GenWcfadsJoin.m \
          $(DEST)/GenWturbPlots.m \
          $(DEST)/GenWmoments.m \
          $(DEST)/GenSfcDiv.m \
          $(DEST)/PlotMultiPanelSample.m \
          $(DEST)/PlotProfs.m \
          $(DEST)/PlotSampleData.m \
          $(DEST)/PlotSampleW250.m \
          $(DEST)/GenCfrac2dHists.m \
          $(DEST)/GenProfPlots.m \
          $(DEST)/GenRain2CloudRatio.m \
          $(DEST)/GenSlope1dPlots.m \
          $(DEST)/GenSlopePlots.m \
          $(DEST)/GenSwupBarCharts.m \
          $(DEST)/GenCfCdBarPlots.m \
          $(DEST)/GenTavgRrateFracAreaDist.m \
          $(DEST)/GenTopSwup.m \
          $(DEST)/MakeSinglePanelFigs.m \
          $(DEST)/PlotDistSet.m \
          $(DEST)/PlotInitialSounding.m \
          $(DEST)/PlotInitialAerosols.m \
          $(DEST)/SetPath.m

CONFIGS = $(DIAG_DIR)/DiagConfig

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)

$(DIAG_DIR)/DiagConfig: DiagConfig
	cp $(<) $(@)

