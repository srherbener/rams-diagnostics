#
# Config
#

DEST = $(HOME)/etc/matlab/common

SCRIPTS = \
          $(DEST)/export_fig.m \
          $(DEST)/figmerge.m \
          $(DEST)/print2array.m \
          $(DEST)/screen2jpeg.m \
          $(DEST)/AdiffEofPlot.m \
          $(DEST)/AzavgDiffEof.m \
          $(DEST)/CompileCounts.m \
          $(DEST)/GenAdiffEofPlots.m \
          $(DEST)/GenCountBins.m \
          $(DEST)/GenDeltapTseries.m \
          $(DEST)/GenHistMeas.m \
          $(DEST)/GenLinePlots.m \
          $(DEST)/GenPcolorPlots.m \
          $(DEST)/GenProfMeas.m \
          $(DEST)/GenSlopes.m \
          $(DEST)/GenTimeAvg.m \
          $(DEST)/GenTsPlots.m \
          $(DEST)/GenVintTseries.m \
          $(DEST)/Hdf5ToBin.m \
          $(DEST)/Hdf5ToNcdf.m \
          $(DEST)/LatToKm.m \
          $(DEST)/LonToKm.m \
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
          $(DEST)/WriteCoordData.m \
          $(DEST)/Xyzt2EofArray.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
