#
# Config
#

DEST = $(HOME)/etc/python/ts_debby

SCRIPTS = \
  $(DEST)/PlotFsFigWind.py \

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

