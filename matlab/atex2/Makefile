#
# Config
#

DEST = $(HOME)/etc/matlab/atex2
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/MarinePrecip/ATEX2

SCRIPTS = \
  $(DEST)/GenAvgFiles.m \
  $(DEST)/SetPath.m

CONFIGS = \
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

