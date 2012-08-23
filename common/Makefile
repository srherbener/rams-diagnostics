#
# Config
#

DEST = $(HOME)/etc/matlab

SCRIPTS = $(DEST)/ReadConfig.m \
          $(DEST)/ReadAzavgVar.m \
          $(DEST)/SmoothFillTseries.m \
          $(DEST)/PlotProfSet.m \
          $(DEST)/PlotTseriesSet.m \
          $(DEST)/Plot2dSet.m \
          $(DEST)/Xyzt2EofArray.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
