dutilsdir = $(prefix)/etc/fortran
dutils_LIBRARIES = libdiag_utils.a

libdiag_utils_a_SOURCES = diag_utils.f90

# Having rhdf5 compile before diag_utils will leave
# the rhdf5_utils.mod file in the rhdf5 (sibling) directory.
AM_FCFLAGS="-I../rhdf5"


DUTILS_MOD_FILE = diag_utils.mod

# Handle the module file
install-exec-local:
	mkdir -p $(dutilsdir)
	cp $(DUTILS_MOD_FILE) $(dutilsdir)

uninstall-local:
	rm -f $(dutilsdir)/$(DUTILS_MOD_FILE)
