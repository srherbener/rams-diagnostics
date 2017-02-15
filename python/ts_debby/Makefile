#
# Config
#

DEST = $(HOME)/etc/python/ts_debby

SCRIPTS = \
  $(DEST)/ConfigTsd.py \
  $(DEST)/GenStormHovmollers.py \
  $(DEST)/FsFigContour.py \
  $(DEST)/PlotFsFigVtFactors.py \
  $(DEST)/PlotFsFigStormXsections.py \
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

