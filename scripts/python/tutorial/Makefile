#
# Config
#

DEST = $(HOME)/etc/python/tutorial

SCRIPTS = \
  $(DEST)/hello.py \
  $(DEST)/test_h5.py \
  $(DEST)/trap.py \
  $(DEST)/simp.py \
  $(DEST)/plot_eg1.py \
  $(DEST)/subplots.py \
  $(DEST)/hdf5.py \
  $(DEST)/hgts_1000mb_map.py \
  $(DEST)/int_examples.py

CONFIGS = 

#
# Targets
#

all: $(DEST) $(SCRIPTS) $(CONFIGS)

$(DEST):
	mkdir -p $(DEST)

$(DEST)/%.py: %.py
	cp $(<) $(@)

$(DIAG_DIR)/%Config: %Config
	cp $(<) $(@)
