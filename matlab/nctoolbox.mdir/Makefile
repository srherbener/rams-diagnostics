#
# Config
#

ML_MAIN = $(HOME)/etc/matlab
DEST = $(ML_MAIN)/nctoolbox

TFILE = /Users/steve/projects/Code/diagnostics/matlab/nctoolbox.mdir/nctoolbox-20130305.tar.gz

#
# Targets
#

all: $(DEST)

$(DEST): $(TFILE)
	mkdir -p $(DEST)
	cd $(DEST); tar -xzvf $(<)
