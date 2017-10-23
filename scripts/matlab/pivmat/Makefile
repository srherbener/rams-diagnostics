#
# Config
#

ML_MAIN = $(HOME)/etc/matlab
DEST = $(ML_MAIN)/pivmat

TFILE = /Users/steve/projects/Code/diagnostics/matlab/pivmat/pivmat.tar.gz

#
# Targets
#

all: $(DEST)

$(DEST): $(TFILE)
	mkdir -p $(DEST)
	cp license.txt $(DEST)
	cd $(DEST); tar -xzvf $(<)
