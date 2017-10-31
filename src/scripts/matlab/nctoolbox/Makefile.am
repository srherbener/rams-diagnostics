NCTB_TARFILE = $(srcdir)/nctoolbox-20130305.tar.gz
NCTB_NAME = nctoolbox
NCTB_INSTALL_DIR = $(prefix)/etc/matlab/$(NCTB_NAME)

EXTRA_DIST = $(NCTB_TARFILE)


# extend main target (all)
all-local: $(NCTB_NAME)

# The touch command is there since the tar command may change the date on
# the $(NCTB_NAME) directory to an older date that causes this target to 
# always appear out of date.
$(NCTB_NAME): $(NCTB_TARFILE)
	mkdir -p $(NCTB_NAME)
	cd $(NCTB_NAME); tar -xzvf ../$(NCTB_TARFILE)
	touch $(NCTB_NAME)

# override install rule since we need a recursive copy, ie use
# cp -r instead of install since install does not do recursive copies
install-exec-local:
	mkdir -p $(NCTB_INSTALL_DIR)
	cp -r $(NCTB_NAME)/* $(NCTB_INSTALL_DIR)

uninstall-local:
	rm -rf $(NCTB_INSTALL_DIR)
