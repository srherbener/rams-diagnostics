#
# Config
#

ML_MAIN = $(HOME)/etc/matlab
DEST = $(ML_MAIN)/m_map

TFILE = /Users/steve/projects/Code/diagnostics/matlab/m_map/m_map1.4.1.tar

#
# Targets
#

all: $(DEST)

$(DEST): $(TFILE)
	cd $(ML_MAIN); tar -xzvf $(<)
