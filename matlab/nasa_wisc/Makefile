#
# Config
#

DEST = $(HOME)/etc/matlab/nasa_wisc
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/NasaWisc

SCRIPTS = \
  $(DEST)/GenDomAvgTseries.m \
  $(DEST)/GenEqMeasTseries.m \
  $(DEST)/Gen2dHovmoller.m \
  $(DEST)/ExtractSampleProfiles.m \
  $(DEST)/ReadTogaCoareSound.m \
  $(DEST)/TogaCoareToRamsSound.m \
  $(DEST)/PlotSampleProfiles.m \
  $(DEST)/SetPath.m

CONFIGS = \
  $(DIAG_DIR)/DiagConfig 

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)

$(DIAG_DIR)/DiagConfig: DiagConfig
	cp $(<) $(@)
