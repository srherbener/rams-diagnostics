# Process this file with autoconf to produce a configure script.
AC_INIT([rams-diagnostics], [1.0.0], [stephen.herbener@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
#AC_CONFIG_SRCDIR([scripts/shell/make_anim_file])
AC_CONFIG_HEADERS([config.h]) 
# Set the default prefix to the user's home directory
AC_PREFIX_DEFAULT([$HOME])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
#AC_CHECK_LIB([c], [_dl_start])
#AC_CHECK_LIB([dl], [dlopen])
#AC_CHECK_LIB([gcc_eh], [_Unwind_Backtrace])
#AC_CHECK_LIB([hdf5], [H5open])
#AC_CHECK_LIB([hdf5_hl], [H5DOread_chunk])
#AC_CHECK_LIB([m], [__ieee754_log])
#AC_CHECK_LIB([pthread], [__pthread_initialize_minimal])
#AC_CHECK_LIB([sz], [szip_check_params])
#AC_CHECK_LIB([z], [compress])

# Checks for header files.
#AC_CHECK_HEADERS([string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([floor gettimeofday mkdir sqrt])

AC_CONFIG_FILES([Makefile
                 lib/Makefile
                 lib/perl/Makefile
                 scripts/Makefile
                 scripts/shell/Makefile
                 scripts/perl/Makefile
                 scripts/idl/Makefile
                 scripts/matlab/Makefile
                 scripts/matlab/atex/Makefile
                 scripts/matlab/atex2/Makefile
                 scripts/matlab/common/Makefile
                 scripts/matlab/ideal_etc/Makefile
                 scripts/matlab/m_map/Makefile
                 scripts/matlab/metpack/Makefile
                 scripts/matlab/nasa_wisc/Makefile
                 scripts/matlab/nctoolbox/Makefile
                 scripts/matlab/obj_anal/Makefile
                 scripts/matlab/pivmat/Makefile
                 scripts/matlab/tc_seed/Makefile
                 scripts/matlab/tempest/Makefile
                 scripts/matlab/ts_debby/Makefile
                 scripts/matlab/vapor/Makefile
                 scripts/python/common/Makefile
                 scripts/python/nasa_wisc/Makefile
                 scripts/python/ts_debby/Makefile
                 scripts/python/tutorial/Makefile
                 scripts/python/Makefile])
AC_OUTPUT

