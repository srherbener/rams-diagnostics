#
# Config
#

DEST = $(HOME)/etc/matlab/rce
DIAG_DIR = /Volumes/tasman1/sherbener/data/RAMS/RCE

SCRIPTS = \
  $(DEST)/GenDomAvgTseries.m \
  $(DEST)/Gen2dHovmoller.m \
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
