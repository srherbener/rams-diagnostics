rhdf5dir = $(prefix)/etc/fortran
rhdf5_LIBRARIES = librhdf5_utils.a

librhdf5_utils_a_SOURCES = rhdf5_f2c.c rhdf5_snames.h rhdf5_utils.f90

RHDF5_MOD_FILE = rhdf5_utils.mod

# Handle the module file
install-exec-local:
	cp $(RHDF5_MOD_FILE) $(rhdf5dir)
