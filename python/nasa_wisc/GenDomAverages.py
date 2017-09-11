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
#    'RCE_S300_SQ',
#    'RCE_S300_SQ_C05',
#    'RCE_S300_SQ_DEBUG',
#    'RCE_S300_SQ_E1',
#    'RCE_S300_SQ_SM',
#    'RCE_1km',
#    'RCE_1km_SM',
#    'RCE_1km_DM',
#    'RCE_1km_DP',
#    'RCE_1km_IPL2',
#    'RCE_1km_SM_IPL2',
#    'RCE_1km_DM_IPL2',
#    'RCE_1km_DP_IPL2',
    'RCE_3km_1mom',
    'RCE_3km_2mom',
    'RCE_3km_2mom_db',
    'RCE_3km_2mom_dm',
    ]
Nsims = len(SimList)

VarList = [
    #  <in_file_prefix> <in_var_name> <out_var_name> <selection_criteria>
    #
    #    <selection_criteria>
    #        "all" - use all points
    #        "ge:<val>" - use points >= to <val>

    [ "lat_flux",  "/lat_flux",  "sfc_lat",  "all" ],
    [ "sens_flux", "/sens_flux", "sfc_sens", "all" ],

    [ "swdn_sfc",   "/rshort",   "sfc_swdn", "all" ],
    [ "swdn_sfc",   "/rshort",   "sfc_swup", "all" ],
    [ "lwdn_sfc",   "/rlong",    "sfc_lwdn", "all" ],
    [ "lwup_sfc",   "/rlongup",  "sfc_lwup", "all" ],

    [ "swdn_toa",  "/swdn",  "top_swdn", "all" ],
    [ "swup_toa",  "/swup",  "top_swup", "all" ],
    [ "lwup_toa",  "/lwup",  "top_lwup", "all" ],

    [ "pcprr",      "/pcprr",         "pcprr",         "all"     ],
    [ "vint_vapor", "/vertint_vapor", "vint_vapor",    "all"     ],
    [ "vint_vapor", "/vertint_vapor", "vint_vapor_01", "ge:0.01" ],
    [ "pi_sfc",     "/pi",            "pi_sfc",        "all" ],

    [ "rtotal", "/rtotal_orig", "total_cond",    "all"     ],
    [ "vapor",  "/vapor",       "vapor",         "all"     ],
    [ "rtotal", "/rtotal_orig", "total_cond_01", "ge:0.01" ],
    [ "vapor",  "/vapor",       "vapor_01",      "ge:0.01" ],
    [ "theta",  "/theta",       "theta",         "all"     ],

    ] 
Nvars = len(VarList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

#InFileTemplate = "HDF5/<SIM>/HDF5/<FPREFIX>-<SIM>-AS-2012-01-01-000000-g1.h5"
InFileTemplate = "SIMDATA/<SIM>/HDF5/<FPREFIX>-<SIM>-LS-2012-01-01-000000-g1.h5"

# Routine to calculate the average of a 2D (horizontal) field
def CalcHorizAvg(Var, Select):
    # Exclude the horizontal edges since RAMS uses those for boundary conditions.
    # Selection process will return a 1D array with all of the selected values
    Var = Var[1:-1, 1:-1]

    if (Select[0] == "all"):
        Var = Var.flatten()
    elif (Select[0] == "ge"):
        Sval = float(Select[1])
        Var = Var[ Var >= Sval ]
        
    return np.mean(Var.astype(np.float64))


# MAIN
for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating horizontal domain averages for simulation: {0:s}".format(Sim))
    print("")

    # Read input variables, calculate horizontal averages and write out results
    for ivar in range(Nvars):
        InFprefix = VarList[ivar][0]
        InVname   = VarList[ivar][1]
        OutVname  = VarList[ivar][2]
        Select    = VarList[ivar][3].split(':')

        # Open the input file
        InFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", InFprefix)
        InFile  = h5py.File(InFname, mode='r')
        print("  Reading: {0:s} ({1:s})".format(InFname, InVname))

        if (OutVname == "sfc_swup"):
            AlbFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "albedt")
            AlbFile  = h5py.File(AlbFname, mode='r')
            AlbVname = "/albedt"
            print("    Reading: {0:s} ({1:s})".format(AlbFname, AlbVname))
        elif (OutVname == "total_cond"):
            RvFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "vapor")
            RvFile  = h5py.File(RvFname, mode='r')
            RvVname = "/vapor"
            print("    Reading: {0:s} ({1:s})".format(RvFname, RvVname))

        # Open the output file
        OutFname = "DIAGS/dom_avg_{0:s}_{1:s}.h5".format(OutVname, Sim)
        OutFile = h5py.File(OutFname, mode='w')

        # Read in coordinates. Try to get Z from a 3D variable.
        print("    Reading: {0:s} ({1:s})".format(InFname, Xname))
        X = InFile[Xname][...]
        Nx = len(X)
        Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
        Xdim.Build(OutFile, X)

        print("    Reading: {0:s} ({1:s})".format(InFname, Yname))
        Y = InFile[Yname][...]
        Ny = len(Y)
        Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
        Ydim.Build(OutFile, Y)

        print("    Reading: {0:s} ({1:s})".format(InFname, Zname))
        Z = InFile[Zname][...]
        Nz = len(Z)
        Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
        Zdim.Build(OutFile, Z)

        print("    Reading: {0:s} ({1:s})".format(InFname, Tname))
        T = InFile[Tname][...]
        Nt = len(T)
        Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
        Tdim.Build(OutFile, T)

        # Read in input variable and form the domain average
        # Do one time step at a time in order to reduce memory usage
        for i in range(Nt):
            # Read in the var and determine if 2D (t,y,x) or 3D (t,z,y,x) field.
            if (OutVname == "sfc_swup"):
                Var = np.squeeze(InFile[InVname][i,...]) * np.squeeze(AlbFile[AlbVname][i,...])
            elif (OutVname == "total_cond"):
                Var = np.squeeze(InFile[InVname][i,...]) - np.squeeze(RvFile[RvVname][i,...])
            else:
                Var = np.squeeze(InFile[InVname][i,...])
            Ndims = Var.ndim

            if (Ndims == 2):
                # Var is (y,x)
                if (i == 0):
                    VarAvg = np.zeros(Nt)
                VarAvg[i] = CalcHorizAvg(Var, Select)
            elif (Ndims == 3):
                # Var is (z,y,x)
                if (i == 0):
                    VarAvg = np.zeros((Nt,Nz))
                for k in range(Nz):
                    VarAvg[i,k] = CalcHorizAvg(np.squeeze(Var[k,:,:]), Select)

            if ((i % 10) == 0):
                print("    Processing time step: {0:d}".format(i))

        print("")

        InFile.close()
        if (OutVname == "sfc_swup"):
            AlbFile.close()
        elif (OutVname == "total_cond"):
            RvFile.close()

        # Write out average
        print("  Writing: {0:s} ({1:s})".format(OutFname, OutVname))
        Dset = h5u.DsetCoards(OutVname, Ndims, VarAvg.shape)
        if (Ndims == 2):
            Dset.Build(OutFile, VarAvg, Tdim)
        elif (Ndims == 3):
            Dset.Build(OutFile, VarAvg, Tdim, Zdim)

        OutFile.close()

        print("")

