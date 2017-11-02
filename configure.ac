# Process this file with autoconf to produce a configure script.
AC_INIT([rams-diagnostics], [1.0.1], [stephen.herbener@gmail.com])

# Place autotools files in subdirectories (help for cleaning up later on)
AC_CONFIG_AUX_DIR([bin])

# Building a foriegn (non-GNU structure) directory structure
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_SRCDIR([src/scripts/shell/make_anim_file])
AC_CONFIG_HEADERS([config.h]) 

# Set the default prefix to the user's home directory
AC_PREFIX_DEFAULT([$HOME])

# Checks for programs.
AC_PROG_CC([gcc])
AC_PROG_FC([pgf90])
AM_PROG_AR
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([c], [close], [], [AC_MSG_ERROR([Unable to find library: libc])])
AC_CHECK_LIB([dl], [dlopen], [], [AC_MSG_ERROR([Unable to find library: libdl])])
AC_CHECK_LIB([sz], [szip_check_params], [], [AC_MSG_ERROR([Unable to find library: libsz])])
AC_CHECK_LIB([z], [compress], [], [AC_MSG_ERROR([Unable to find library: libz])])
AC_CHECK_LIB([m], [cos], [], [AC_MSG_ERROR([Unable to find library: libm])])
AC_CHECK_LIB([pthread], [pthread_create], [], [AC_MSG_ERROR([Unable to find library: libpthread])])
AC_CHECK_LIB([hdf5], [H5open], [], [AC_MSG_ERROR([Unable to find library: libhdf5])])
AC_CHECK_LIB([hdf5_hl], [H5LTmake_dataset_float], [], [AC_MSG_ERROR([Unable to find library: libhdf5_hl])])

# Checks for header files.
AC_CHECK_HEADERS([string.h sys/time.h], [], [AC_MSG_ERROR([Unable to find include files: string.h, ssy/tim.h])])
AC_CHECK_HEADERS([hdf5.h], [], [AC_MSG_ERROR([Unable to find include file: hdf5.h])])

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([floor gettimeofday mkdir sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/scripts/Makefile
                 src/scripts/perl/Makefile
                 src/scripts/shell/Makefile
                 src/scripts/idl/Makefile
                 src/scripts/matlab/Makefile
                 src/scripts/matlab/atex/Makefile
                 src/scripts/matlab/atex2/Makefile
                 src/scripts/matlab/common/Makefile
                 src/scripts/matlab/ideal_etc/Makefile
                 src/scripts/matlab/m_map/Makefile
                 src/scripts/matlab/metpack/Makefile
                 src/scripts/matlab/nasa_wisc/Makefile
                 src/scripts/matlab/nctoolbox/Makefile
                 src/scripts/matlab/obj_anal/Makefile
                 src/scripts/matlab/pivmat/Makefile
                 src/scripts/matlab/tc_seed/Makefile
                 src/scripts/matlab/tempest/Makefile
                 src/scripts/matlab/ts_debby/Makefile
                 src/scripts/matlab/vapor/Makefile
                 src/scripts/python/common/Makefile
                 src/scripts/python/nasa_wisc/Makefile
                 src/scripts/python/ts_debby/Makefile
                 src/scripts/python/tutorial/Makefile
                 src/scripts/python/Makefile
                 src/libraries/Makefile
                 src/libraries/getoptions/Makefile
                 src/libraries/rhdf5/Makefile
                 src/libraries/diag_utils/Makefile
                 src/programs/Makefile
                 src/programs/azavg/Makefile
                 src/programs/create_bins/Makefile
                 src/programs/diag_filter/Makefile
                 src/programs/find_tc_center/Makefile
                 src/programs/gen_adv/Makefile
                 src/programs/gen_flux/Makefile
                 src/programs/gen_moments/Makefile
                 src/programs/gen_tc_tilt/Makefile
                 src/programs/gen_u_v/Makefile
                 src/programs/gen_vint_terms/Makefile
                 src/programs/gen_vt_vr/Makefile
                 src/programs/gen_xsection/Makefile
                 src/programs/hdata_op/Makefile
                 src/programs/join_hdata/Makefile
                 src/programs/sig_proc/Makefile
                 src/programs/sub_vortex/Makefile
                 src/programs/tsavg/Makefile])


AC_OUTPUT

