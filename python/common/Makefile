#
# Config
#

DEST = $(HOME)/etc/python/common

SCRIPTS = \
  $(DEST)/Hdf5Utils.py \
  $(DEST)/PlotUtils.py \
  $(DEST)/PlotTcXsection.py

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

