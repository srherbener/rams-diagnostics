#
# Config
#

DEST = $(HOME)/etc/matlab/metpack

SCRIPTS = \
          $(DEST)/getcape.m \
          $(DEST)/lcl.m \
          $(DEST)/skew_sounding.m \
          $(DEST)/tskew.m \
          $(DEST)/windbarb.m

#
# Targets
#

all: $(DEST) $(SCRIPTS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.m: %.m
	cp $(<) $(@)
