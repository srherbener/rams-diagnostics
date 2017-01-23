#
# Config
#

DEST = $(HOME)/etc/python/tutorial

SCRIPTS = \
  $(DEST)/hello.py \
  $(DEST)/trap.py \
  $(DEST)/simp.py \
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
