#
# Config
#

DEST = $(HOME)/etc/python/nasa_wisc

SCRIPTS = \
  $(DEST)/ConfigRce.py \
  $(DEST)/ExtractRceLiteVars.py \
  $(DEST)/GenRamsDomAvgs.py \
  $(DEST)/GenDomAverages.py \
  $(DEST)/GenEqMeas.py \
  $(DEST)/GenProfiles.py \
  $(DEST)/PlotProfiles.py \
  $(DEST)/PlotTimeSeries.py


#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

