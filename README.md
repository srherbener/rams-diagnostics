# rams-diagnostics

This repository contains Fortran programs and scripts in various languages (bash, MATLAB, Python, Perl).

This software is used for analyzing the output of RAMS (Regional Atmospheric Modeling System), or REVU which is
the dataset extractor for RAMS. For example, the FORTRAN program azavg will do azimuthal averaging of atmospheric
quantities in simulated tropical cyclones.

The programs use the C interface to HDF5 and thus need to be told the paths to the HDF5 include files and
and libraries. This is accomplished by passing in CPPFLAGS and LDFLAGS settings when running configure.

For example, let's say you have installed the HDF5 package in the path /opt/hdf5. Then run configure
as follows:

$ ./configure LDFLAGS="-L/opt/hdf5/lib" CPPFLAGS="-I/opt/path/include"

The default prefix in this package has been set to $HOME. Use --prefix=new_path to override.

