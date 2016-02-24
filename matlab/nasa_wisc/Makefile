#
# Config
#

DEST = $(HOME)/etc/matlab/nasa_wisc
DIAG_DIR = $(HOME)/projects/NasaWisc/Work

SCRIPTS = \
  $(DEST)/GenDomAvgTseries.m \
  $(DEST)/GenEqMeasTseries.m \
  $(DEST)/Gen2dHovmoller.m \
  $(DEST)/GenMeanHydroParams.m \
  $(DEST)/GenRceHistMeas.m \
  $(DEST)/ExtractSampleProfiles.m \
  $(DEST)/ReadTogaCoareSound.m \
  $(DEST)/TogaCoareToRamsSound.m \
  $(DEST)/PlotSampleProfiles.m \
  $(DEST)/PlotEqMeas.m \
  $(DEST)/SetPath.m

CONFIGS = \
  $(DIAG_DIR)/DiagConfig \
  $(DIAG_DIR)/FigLatHeatProfileConfig \
  $(DIAG_DIR)/FigureConfig \
  $(DIAG_DIR)/FigReport062415_Config

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)

$(DIAG_DIR)/%Config: %Config
	cp $(<) $(@)
