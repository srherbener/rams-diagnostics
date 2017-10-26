# Need this for the substitutions in rhdf5_snames.h
AM_CFLAGS = -DPC_LINUX1

rhdf5dir = $(prefix)/etc/fortran
rhdf5_LIBRARIES = librhdf5_utils.a

librhdf5_utils_a_SOURCES = rhdf5_f2c.c rhdf5_snames.h rhdf5_utils.f90

RHDF5_MOD_FILE = rhdf5_utils.mod

# Handle the module file
install-exec-local:
	mkdir -p $(rhdf5dir)
	cp $(RHDF5_MOD_FILE) $(rhdf5dir)

uninstall-local:
	rm -f $(rhdf5dir)/$(RHDF5_MOD_FILE)
