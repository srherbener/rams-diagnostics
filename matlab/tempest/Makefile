#
# Config
#

DEST = $(HOME)/etc/matlab/tempest
DIAG_DIR = /Volumes/data/RAMS/NasaTempest

SCRIPTS = $(DEST)/GenTotPcp.m \
          $(DEST)/SetPath.m \
          $(DEST)/Sounding.m \
          $(DEST)/Sounding2.m

CONFIGS = $(DIAG_DIR)/DiagConfig

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

