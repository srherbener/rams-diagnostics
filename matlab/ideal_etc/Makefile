#
# Config
#

DEST = $(HOME)/etc/matlab/ideal_etc
#DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/MarinePrecip/ATEX

SCRIPTS = $(DEST)/GenInitFields.m \
          $(DEST)/SetPath.m

#CONFIGS = $(DIAG_DIR)/DiagConfig
CONFIGS =

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)

#$(DIAG_DIR)/DiagConfig: DiagConfig
#	cp $(<) $(@)

