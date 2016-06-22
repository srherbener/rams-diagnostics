#
# Config
#

DEST = $(HOME)/etc/matlab/common

SCRIPTS = \
          $(DEST)/export_fig.m \
          $(DEST)/freezeColors.m \
          $(DEST)/unfreezeColors.m \
          $(DEST)/cbfreeze.m \
          $(DEST)/cbhandle.m \
          $(DEST)/test_freezeColors.m \
          $(DEST)/ReduceData.m \
          $(DEST)/SetPlotAxes.m \
          $(DEST)/DrawPlotData.m \
          $(DEST)/GenFigures.m \
          $(DEST)/GenAxesSpecs.m \
          $(DEST)/GenDataSpecs.m \
          $(DEST)/GenLinePlot.m \
          $(DEST)/LeftJustTitle.m \
          $(DEST)/AssociateStructs.m \
          $(DEST)/Line2Fields.m \
          $(DEST)/ParseTitle.m \
          $(DEST)/ReadFigureConfig.m \
          $(DEST)/ReadPrepConfigFile.m \
          $(DEST)/WriteStringAttribute.m \
          $(DEST)/CreateDimensionsXyzt.m \
          $(DEST)/AttachDimensionsXyzt.m \
          $(DEST)/NotateDimensionsXyzt.m \
          $(DEST)/NotateVariableXyzt.m \
          $(DEST)/figmerge.m \
          $(DEST)/gridLegend.m \
          $(DEST)/columnlegend.m \
          $(DEST)/inpaintn.m \
          $(DEST)/print2array.m \
          $(DEST)/str2rgb.m \
          $(DEST)/screen2jpeg.m \
          $(DEST)/radar_cmap.m \
          $(DEST)/AdiffEofPlot.m \
          $(DEST)/AzavgDiffEof.m \
          $(DEST)/CombineBins.m \
          $(DEST)/CompileCounts.m \
          $(DEST)/CountsToAvg.m \
          $(DEST)/Gen2dHistData.m \
          $(DEST)/GenAdiffEofPlots.m \
          $(DEST)/GenCountBins1d.m \
          $(DEST)/GenCountBins.m \
          $(DEST)/GenDeltapTseries.m \
          $(DEST)/GenHistMeas.m \
          $(DEST)/GenLinePlots.m \
          $(DEST)/GenMoments.m \
          $(DEST)/GenContourPlots.m \
          $(DEST)/GenPcolorPlots.m \
          $(DEST)/GenProfMeas.m \
          $(DEST)/GenSlopes1d.m \
          $(DEST)/GenSlopes.m \
          $(DEST)/GenSoundingStruct.m \
          $(DEST)/GenTimeAvg.m \
          $(DEST)/GenTsPlots.m \
          $(DEST)/GenVintTseries.m \
          $(DEST)/Hdf5ToBin.m \
          $(DEST)/Hdf5ToNcdf.m \
          $(DEST)/LatToKm.m \
          $(DEST)/LonToKm.m \
          $(DEST)/MakeBlankPlot.m \
          $(DEST)/ReadConfig.m \
          $(DEST)/ReadAzavgVar.m \
          $(DEST)/ReadSelectXyzt.m \
          $(DEST)/ReadSelectHdf5.m \
          $(DEST)/ReduceHists.m \
          $(DEST)/ReduceHist1D.m \
          $(DEST)/RunAzavgDiffEof.m \
          $(DEST)/SmoothFillTseries.m \
          $(DEST)/SoundingToProfile.m \
          $(DEST)/PlotProfSet.m \
          $(DEST)/PlotTseriesSet.m \
          $(DEST)/Plot2dSet.m \
          $(DEST)/PlotBarSet.m \
          $(DEST)/PlotContour.m \
          $(DEST)/TileFigs.m \
          $(DEST)/TimeToString.m \
          $(DEST)/WriteCoordData.m \
          $(DEST)/WriteDpFile.m \
          $(DEST)/ReadDpFile.m \
          $(DEST)/Xyzt2EofArray.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
