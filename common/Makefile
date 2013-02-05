#
# Config
#

DEST = $(HOME)/etc/matlab

SCRIPTS = \
          $(DEST)/figmerge.m \
          $(DEST)/AdiffEofPlot.m \
          $(DEST)/AzavgDiffEof.m \
          $(DEST)/GenAdiffEofPlots.m \
          $(DEST)/GenDeltapTseries.m \
          $(DEST)/GenHistMeas.m \
          $(DEST)/GenLinePlots.m \
          $(DEST)/GenProfMeas.m \
          $(DEST)/GenTimeAvg.m \
          $(DEST)/GenTsPlots.m \
          $(DEST)/GenVintTseries.m \
          $(DEST)/ReadConfig.m \
          $(DEST)/ReadAzavgVar.m \
          $(DEST)/ReadSelectXyzt.m \
          $(DEST)/ReduceHists.m \
          $(DEST)/RunAzavgDiffEof.m \
          $(DEST)/SmoothFillTseries.m \
          $(DEST)/PlotProfSet.m \
          $(DEST)/PlotTseriesSet.m \
          $(DEST)/Plot2dSet.m \
          $(DEST)/TileFigs.m \
          $(DEST)/TimeToString.m \
          $(DEST)/Xyzt2EofArray.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
