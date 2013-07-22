#
# Config
#

DEST = $(HOME)/etc/perl

SCRIPTS = $(DEST)/DiagUtils.pm

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.pl: %.pl
	cp $(<) $(@)

$(DEST)/%.pm: %.pm
	cp $(<) $(@)
