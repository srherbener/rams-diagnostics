PM_TARFILE = $(srcdir)/pivmat.tar.gz
PM_NAME = pivmat
PM_INSTALL_DIR = $(prefix)/etc/matlab/$(PM_NAME)

EXTRA_DIST = $(PM_TARFILE)


# extend main target (all)
all-local: $(PM_NAME)

# The touch command is there since the tar command may change the date on
# the $(PM_NAME) directory to an older date that causes this target to 
# always appear out of date.
$(PM_NAME): $(PM_TARFILE)
	mkdir -p $(PM_NAME)
	cd $(PM_NAME); tar -xzvf ../$(PM_TARFILE)
	touch $(PM_NAME)

# override install rule since we need a recursive copy, ie use
# cp -r instead of install since install does not do recursive copies
install-exec-local:
	mkdir -p $(PM_INSTALL_DIR)
	cp -r $(PM_NAME)/* $(PM_INSTALL_DIR)

uninstall-local:
	rm -rf $(PM_INSTALL_DIR)
