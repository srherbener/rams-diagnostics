#
# Config
#

DEST = $(HOME)/etc/matlab/tc_seed

SCRIPTS = $(DEST)/EofAnalysis.m \
          $(DEST)/SetPath.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
