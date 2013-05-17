#
# Config
#

DEST = $(HOME)/etc/matlab/metpack

SCRIPTS = \
          $(DEST)/getcape.m \
          $(DEST)/lcl.m \
          $(DEST)/skew_sounding.m \
          $(DEST)/soundings_RF6.mat \
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

$(DEST)/soundings_RF6.mat: soundings_RF6.mat
	cp $(<) $(@)
