#
# Config
#

DEST = $(HOME)/etc/python/ts_debby

SCRIPTS = \
  $(DEST)/ConfigTsd.py \
  $(DEST)/PlotFsFigVtFactors.py \
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

