#
# Config
#

ML_MAIN = $(HOME)/etc/matlab
DEST = $(ML_MAIN)/m_map

TFILE = $(PWD)/m_map1.4.tar

#
# Targets
#

all: $(DEST)

$(DEST): $(TFILE)
	cd $(ML_MAIN); tar -xzvf $(<)

$(DEST)/%.m: %.m
	cp $(<) $(@)
