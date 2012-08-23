#
# Config
#

DEST = $(HOME)/etc/matlab/atex

SCRIPTS = $(DEST)/GenProfPlots.m \
          $(DEST)/GenSwPlots.m \
          $(DEST)/PlotProfSet.m \
          $(DEST)/SetPath.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
