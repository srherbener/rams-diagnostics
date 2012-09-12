#
# Config
#

DEST = $(HOME)/etc/matlab/atex

SCRIPTS = $(DEST)/DomainTimeAvgProfile.m \
          $(DEST)/ExtractSample.m \
          $(DEST)/GenDomAvg.m \
          $(DEST)/GenDomTimAvgProfs.m \
          $(DEST)/GenHist2d.m \
          $(DEST)/GenProfPlots.m \
          $(DEST)/GenRrateHist.m \
          $(DEST)/GenSwPlots.m \
          $(DEST)/GenTimeSeries.m \
          $(DEST)/MergeHists.m \
          $(DEST)/PlotHist2d.m \
          $(DEST)/PlotHistSet.m \
          $(DEST)/PlotHists.m \
          $(DEST)/PlotMultiPanelSample.m \
          $(DEST)/PlotRrateHist.m \
          $(DEST)/PlotSampleData.m \
          $(DEST)/PlotTimeSeries.m \
          $(DEST)/RatioAvg2d.m \
          $(DEST)/SetPath.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
