#
# Config
#

DEST = $(HOME)/etc/matlab/atex2
DIAG_DIR = $(HOME)/projects/MarinePrecip/Work/ATEX2

SCRIPTS = \
  $(DEST)/GenAvgFiles.m \
  $(DEST)/GenBarGraphFiles.m \
  $(DEST)/GenBuoyancyFlux.m \
  $(DEST)/GenPdfFiles.m \
  $(DEST)/GenMomentFiles.m \
  $(DEST)/SetPath.m

CONFIGS = \
  $(DIAG_DIR)/FigMicrophysicsS293_Config \
  $(DIAG_DIR)/FigMicrophysicsS298_Config \
  $(DIAG_DIR)/FigProfilesS293_Config \
  $(DIAG_DIR)/FigProfilesS298_Config \
  $(DIAG_DIR)/FigBarGraphs_Config \
  $(DIAG_DIR)/FigBlStats_Config \
  $(DIAG_DIR)/FigTestInitCcn_Config \
  $(DIAG_DIR)/FigTestTimeSeries_Config

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

