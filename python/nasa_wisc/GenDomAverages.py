#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigRce as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()

SimList = [
    'RCE_S300_SQ',
    'RCE_S300_SQ_C05',
    ]
Nsims = len(SimList)

VarList = [
    #  <in_file_prefix> <in_var_name> <out_var_name>
    [ "lat_flux",  "/lat_flux",  "/avg_sfc_lat"  ],
    [ "sens_flux", "/sens_flux", "/avg_sfc_sens" ],

    [ "rshort",   "/rshort",   "/avg_sfc_swdn" ],
    [ "rshortup", "/rshortup", "/avg_sfc_swup" ],
    [ "rlong",    "/rlong",    "/avg_sfc_lwdn" ],
    [ "rlongup",  "/rlongup",  "/avg_sfc_lwup" ],

    [ "albedt", "/albedt", "/avg_sfc_alb"  ],

    [ "swdn_toa",  "/swdntop",  "/avg_top_swdn" ],
    [ "swup_toa",  "/swuptop",  "/avg_top_swup" ],
    [ "lwup_toa",  "/olr",      "/avg_top_lwup" ],

    [ "pcprate",  "/pcprate", "/avg_pcprate" ],

    [ "vint_cond",  "/vertint_cond",  "/avg_vint_cond"  ],
    [ "vint_vapor", "/vertint_vapor", "/avg_vint_vapor" ],
    ] 
Nvars = len(VarList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

InFileTemplate = "HDF5/<SIM>/HDF5/<FPREFIX>-<SIM>-AS-2012-01-01-000000-g1.h5"

for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating horizontal domain averages for simulation: {0:s}".format(Sim))
    print("")

    # Open the output file
    OutFname = "DIAGS/dom_avgerage_{0:s}.h5".format(Sim)
    OutFile = h5py.File(OutFname, mode='w')

    # Read input variables, calculate horizontal averages and write out results
    for ivar in range(Nvars):
        InFprefix = VarList[ivar][0]
        InVname   = VarList[ivar][1]
        OutVname  = VarList[ivar][2]

        # Open the input file
        InFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", InFprefix)
        InFile  = h5py.File(InFname, mode='r')

        # If first variable, read in the coordinates and copy to the output file
        if (ivar == 0):
            print("  Reading: {0:s} ({1:s})".format(InFname, Xname))
            print("  Reading: {0:s} ({1:s})".format(InFname, Yname))
            print("  Reading: {0:s} ({1:s})".format(InFname, Zname))
            print("  Reading: {0:s} ({1:s})".format(InFname, Tname))
            print("")

            X = InFile[Xname][...]
            Y = InFile[Yname][...]
            Z = InFile[Zname][...]
            T = InFile[Tname][...]

            Nx = len(X)
            Ny = len(Y)
            Nz = len(Z)
            Nt = len(T)

            # Write out the coordinate data, mark these as dimensions
            Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
            Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
            Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
            Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)

            Xdim.Build(OutFile, X)
            Ydim.Build(OutFile, Y)
            Zdim.Build(OutFile, Z)
            Tdim.Build(OutFile, T)

        # Read in input variable and form the domain average
        # Do one time step at a time in order to reduce memory usage
        print("  Reading: {0:s} ({1:s})".format(InFname, InVname))
        VarAvg  = np.zeros(Nt)
        for i in range(Nt):
            # Exclude the horizontal edges since RAMS uses those for boundary conditions.
            Var = np.squeeze(InFile[InVname][i,1:-1,1:-1])
            VarAvg[i] = np.mean(Var.astype(np.float64))

            if ((i % 10) == 0):
                print("    Processing time step: {0:d}".format(i))

        print("")

        InFile.close()

        # Write out average
        print("  Writing: {0:s} ({1:s})".format(OutFname, OutVname))
        Dset = h5u.DsetCoards(OutVname, 1, [ Nt ])
        Dset.Build(OutFile, VarAvg, Tdim)

        print("")

    OutFile.close()
