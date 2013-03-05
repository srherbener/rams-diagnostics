#
# Config
#

DEST = $(HOME)/etc/matlab/atex
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/MarinePrecip/ATEX

SCRIPTS = $(DEST)/ExtractSample.m \
          $(DEST)/GenCloudStructPlots.m \
          $(DEST)/GenDistPlots.m \
          $(DEST)/PlotMultiPanelSample.m \
          $(DEST)/PlotSampleData.m \
          $(DEST)/GenProfPlots.m \
          $(DEST)/GenRain2CloudRatio.m \
          $(DEST)/GenSlopePlots.m \
          $(DEST)/GenSwupBarCharts.m \
          $(DEST)/GenTavgRrateFracAreaDist.m \
          $(DEST)/GenTopSwup.m \
          $(DEST)/MakeSinglePanelFigs.m \
          $(DEST)/PlotDistSet.m \
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

