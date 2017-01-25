#
# Config
#

DEST = $(HOME)/etc/python/common

SCRIPTS = \
  $(DEST)/PlotUtils.py \

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

