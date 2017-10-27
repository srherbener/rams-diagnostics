MMAP_TARFILE = $(srcdir)/m_map1.4.1.tar
MMAP_NAME = m_map
MMAP_INSTALL_DIR = $(prefix)/etc/matlab/$(MMAP_NAME)

EXTRA_DIST = $(MMAP_TARFILE)


# extend main target (all)
all-local: $(MMAP_NAME)

# The touch command is there since the tar command may change the date on
# the $(MMAP_NAME) directory to an older date that causes this target to 
# always appear out of date.
$(MMAP_NAME): $(MMAP_TARFILE)
	mkdir -p $(MMAP_NAME)
	tar -xzvf $(MMAP_TARFILE)
	touch $(MMAP_NAME)

# override install rule since we need a recursive copy, ie use
# cp -r instead of install since install does not do recursive copies
install-exec-local:
	mkdir -p $(MMAP_INSTALL_DIR)
	cp -r $(MMAP_NAME)/* $(MMAP_INSTALL_DIR)

uninstall-local:
	rm -rf $(MMAP_INSTALL_DIR)
