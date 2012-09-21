#
# Config
#

DEST = $(HOME)/etc/matlab/atex

SCRIPTS = $(DEST)/ExtractSample.m \
          $(DEST)/GenCloudStructPlots.m \
          $(DEST)/GenDistPlots.m \
          $(DEST)/PlotMultiPanelSample.m \
          $(DEST)/PlotSampleData.m \
          $(DEST)/GenProfPlots.m \
          $(DEST)/GenRain2CloudRatio.m \
          $(DEST)/GenSwupBarCharts.m \
          $(DEST)/GenTavgRrateFracAreaDist.m \
          $(DEST)/GenTopSwup.m \
          $(DEST)/PlotDistSet.m \
          $(DEST)/SetPath.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
