#
# Config
#

DEST = $(HOME)/etc/matlab/atex
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/MarinePrecip/ATEX

SCRIPTS = $(DEST)/ExtractSample.m \
          $(DEST)/GenAvgCloudDepth.m \
          $(DEST)/GenAlbedo.m \
          $(DEST)/GenBLcouplingData.m \
          $(DEST)/GenBuoyancyFlux.m \
          $(DEST)/GenBlStats.m \
          $(DEST)/GenCloudStructPlots.m \
          $(DEST)/GenDistPlots.m \
          $(DEST)/GenDropSizeDists.m \
          $(DEST)/GenHfracCloud.m \
          $(DEST)/GenBarGraphData.m \
          $(DEST)/GenMomentData.m \
          $(DEST)/GenProfileData.m \
          $(DEST)/GenHistData.m \
          $(DEST)/GenPdfPlots.m \
          $(DEST)/GenTurbData.m \
          $(DEST)/GenWcfads.m \
          $(DEST)/GenWcfadsJoin.m \
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

