#
# Config
#

DEST = $(HOME)/etc/python/common

SCRIPTS = \
  $(DEST)/Hdf5Utils.py \
  $(DEST)/GetDiskUsage.py \
  $(DEST)/PlotUtils.py

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

