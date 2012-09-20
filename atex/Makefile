#
# Config
#

DEST = $(HOME)/etc/matlab/atex

SCRIPTS = $(DEST)/GenCloudStructPlots.m \
          $(DEST)/GenDistPlots.m \
          $(DEST)/GenProfPlots.m \
          $(DEST)/GenRain2CloudRatio.m \
          $(DEST)/GenTavgRrateFracAreaDist.m \
          $(DEST)/GenTopSwup.m \
          $(DEST)/PlotDistSet.m \
          $(DEST)/SetPath.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
