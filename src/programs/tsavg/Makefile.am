bin_PROGRAMS = tsavg

tsavg_SOURCES = tsavg.f90

# Having rhdf5 (libraries) compile before this program will place
# the library and module files in the build directory.
MY_LIB_ROOT = ../../libraries
AM_FCFLAGS = -I$(MY_LIB_ROOT)/getoptions -I$(MY_LIB_ROOT)/rhdf5 -I$(MY_LIB_ROOT)/diag_utils
AM_LDFLAGS = -L$(MY_LIB_ROOT)/getoptions -L$(MY_LIB_ROOT)/rhdf5 -L$(MY_LIB_ROOT)/diag_utils
tsavg_LDADD = -lgetoptions -ldiag_utils -lrhdf5_utils
