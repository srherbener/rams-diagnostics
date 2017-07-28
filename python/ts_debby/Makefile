#
# Config
#

DEST = $(HOME)/etc/python/ts_debby

SCRIPTS = \
  $(DEST)/ConfigTsd.py \
  $(DEST)/GenReduceData.py \
  $(DEST)/GenStormHovmollers.py \
  $(DEST)/GenTrackAverages.py \
  $(DEST)/GenPtrackTempGrads.py \
  $(DEST)/GenPtrackMomentum.py \
  $(DEST)/GenFieldNoVortex.py \
  $(DEST)/GenVentIndex.py \
  $(DEST)/FsFigContour.py \
  $(DEST)/PlotFsFigStormHovmollers.py \
  $(DEST)/PlotFsFigStormXsections.py \
  $(DEST)/PlotFsFigTrackXsections.py \
  $(DEST)/PlotFsFigPtrackTempGrad.py \
  $(DEST)/PlotFsFigPlanViews.py \
  $(DEST)/PlotFsFigWind.py \
  $(DEST)/PlotFsFigWindTseries.py

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

