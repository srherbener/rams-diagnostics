#
# Config
#

DEST = $(HOME)/etc/python/ts_debby

SCRIPTS = \
  $(DEST)/ConfigTsd.py \
  $(DEST)/GenStormHovmollers.py \
  $(DEST)/GenPtrackAverages.py \
  $(DEST)/GenPtrackTempGrads.py \
  $(DEST)/GenPtrackMomentum.py \
  $(DEST)/GenWindNoVortex.py \
  $(DEST)/FsFigContour.py \
  $(DEST)/PlotFsFigVtFactors.py \
  $(DEST)/PlotFsFigStormHovmollers.py \
  $(DEST)/PlotFsFigStormXsections.py \
  $(DEST)/PlotFsFigTrackXsections.py \
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

