#
# Config
#

DEST = $(HOME)/etc/matlab

SCRIPTS = $(DEST)/ReadConfig.m \
          $(DEST)/PlotTseriesSet.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
