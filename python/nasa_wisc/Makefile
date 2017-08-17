#
# Config
#

DEST = $(HOME)/etc/python/nasa_wisc

SCRIPTS = \
  $(DEST)/ConfigRce.py \
  $(DEST)/ExtractRceLiteVars.py \
  $(DEST)/GenRamsDomAvgs.py \
  $(DEST)/GenDomAverages.py \
  $(DEST)/GenEqMeas.py


#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

