#
# Config
#

DEST = $(HOME)/etc/python/nasa_wisc

SCRIPTS = \
  $(DEST)/ConfigRce.py \
  $(DEST)/ExtractRceLiteVars.py \
  $(DEST)/GenRamsDomAvgs.py \
  $(DEST)/GenDomAverages.py \
  $(DEST)/GenAnimFrames.py \
  $(DEST)/GenEqMeas.py \
  $(DEST)/GenProfiles.py \
  $(DEST)/GenPowerSpectra.py \
  $(DEST)/PlotProfiles.py \
  $(DEST)/PlotPowerSpectra.py \
  $(DEST)/PlotTimeSeries.py


#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

