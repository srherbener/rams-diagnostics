#
# Config
#

DEST = $(HOME)/etc/matlab/tc_seed

SCRIPTS = $(DEST)/GenTsPlots.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
