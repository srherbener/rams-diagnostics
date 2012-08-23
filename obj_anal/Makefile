#
# Config
#

DEST = $(HOME)/etc/matlab/obj_anal

SCRIPTS = $(DEST)/CompositeEiTemp.m \
          $(DEST)/CovNan.m \
          $(DEST)/Create2dMap.m \
          $(DEST)/CreateCyclicIndices.m \
          $(DEST)/CreateCyclicPattern.m \
          $(DEST)/DiffOfMeansTstat.m \
          $(DEST)/DotProdNan.m \
          $(DEST)/EigenSpectrum.m \
          $(DEST)/EofNan.m \
          $(DEST)/GenCoherence.m \
          $(DEST)/GenCrossSpectrum.m \
          $(DEST)/GenPlotEigSpectrum.m \
          $(DEST)/GenPowerSpectrum.m \
          $(DEST)/GenRnoiseSpectrum.m \
          $(DEST)/GenSampleStats.m \
          $(DEST)/GenSpectrum.m \
          $(DEST)/GenSplitCrossSpectrum.m \
          $(DEST)/GenSplitSpectrum.m \
          $(DEST)/HannWindow.m \
          $(DEST)/LagRegCor.m \
          $(DEST)/LagRegFit.m \
          $(DEST)/LowPassButter.m \
          $(DEST)/MatMultNan.m \
          $(DEST)/MeanNan.m \
          $(DEST)/NormCols.m \
          $(DEST)/Plot2dEofPc.m \
          $(DEST)/Plot2dMap.m \
          $(DEST)/PlotPowerSpectrum.m \
          $(DEST)/PlotWorldDataMap.m \
          $(DEST)/PrintSigPeaks.m \
          $(DEST)/PullOutCyclicData.m \
          $(DEST)/ReadEnsoFile.m \
          $(DEST)/ReadGlobalTempFile.m \
          $(DEST)/RedNoise.m \
          $(DEST)/RedNoiseSpectrum.m \
          $(DEST)/RegCorCoeffs.m \
          $(DEST)/RegFit.m \
          $(DEST)/RemoveMeanNan.m \
          $(DEST)/RunFft.m \
          $(DEST)/RunningMean.m \
          $(DEST)/SetAreaToOne.m \
          $(DEST)/StandardizeData.m \
          $(DEST)/SweepLagCor.m \
          $(DEST)/VarNan.m \
          $(DEST)/WeightMean.m \
          $(DEST)/myFact.m \
          $(DEST)/redblue.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
