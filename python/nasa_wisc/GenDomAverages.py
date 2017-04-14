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
    'RCE_S300_SQ_DEBUG',
    'RCE_S300_SQ_E1',
    'RCE_S300_SQ_SM',
    ]
Nsims = len(SimList)

VarList = [
    #  <in_file_prefix> <in_var_name> <out_var_name> <selection_criteria>
    #
    #    <selection_criteria>
    #        "all" - use all points
    #        "ge:<val>" - use points >= to <val>

    [ "lat_flux",  "/lat_flux",  "/avg_sfc_lat",  "all" ],
    [ "sens_flux", "/sens_flux", "/avg_sfc_sens", "all" ],

    [ "rshort",   "/rshort",   "/avg_sfc_swdn", "all" ],
    [ "rshortup", "/rshortup", "/avg_sfc_swup", "all" ],
    [ "rlong",    "/rlong",    "/avg_sfc_lwdn", "all" ],
    [ "rlongup",  "/rlongup",  "/avg_sfc_lwup", "all" ],

    [ "albedt", "/albedt", "/avg_sfc_alb", "all" ],

    [ "swdn_tod",  "/swdn",  "/avg_top_swdn", "all" ],
    [ "swup_tod",  "/swup",  "/avg_top_swup", "all" ],
    [ "lwup_tod",  "/lwup",  "/avg_top_lwup", "all" ],

    [ "pcprate",  "/pcprate", "/avg_pcprate", "all" ],

    [ "vint_cond",  "/vertint_cond",  "/avg_vint_cond",  "all" ],
    [ "vint_vapor", "/vertint_vapor", "/avg_vint_vapor", "all" ],
    ] 
Nvars = len(VarList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

InFileTemplate = "HDF5/<SIM>/HDF5/<FPREFIX>-<SIM>-AS-2012-01-01-000000-g1.h5"

# Routine to calculate the average of a 2D (horizontal) field
def CalcHorizAvg(Var, Select):
    # Exclude the horizontal edges since RAMS uses those for boundary conditions.
    # Selection process will return a 1D array with all of the selected values
    Var = Var[1:-1, 1:-1]

    if (Select[0] == "all"):
        Var = Var.flatten()
    if (Select[0] == "ge"):
        Sval = float(Select[1])
        Var = Var[ Var >= Sval ]
        
    return np.mean(Var.astype(np.float64))


# MAIN
for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating horizontal domain averages for simulation: {0:s}".format(Sim))
    print("")

    # Open the output file
    OutFname = "DIAGS/dom_averages_{0:s}.h5".format(Sim)
    OutFile = h5py.File(OutFname, mode='w')

    NoX = True
    NoY = True
    NoZ = True
    NoT = True

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

        # Read in coordinates. Try to get Z from a 3D variable.
        if (NoX):
            print("    Reading: {0:s} ({1:s})".format(InFname, Xname))
            X = InFile[Xname][...]
            Nx = len(X)
            Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
            Xdim.Build(OutFile, X)
            NoX = False

        if (NoY):
            print("    Reading: {0:s} ({1:s})".format(InFname, Yname))
            Y = InFile[Yname][...]
            Ny = len(Y)
            Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
            Ydim.Build(OutFile, Y)
            NoY = False

        if (NoZ):
            print("    Reading: {0:s} ({1:s})".format(InFname, Zname))
            Z = InFile[Zname][...]
            Nz = len(Z)
            # Keep grabbing Z coordinates until one cooresponding to a 3D field
            # is obtained. If a 3D Z coordinate vector is found, write it
            # out. If one is not found, GRADS is okay with the z coordinates
            # missing, plus all the variables are 2D and the z coordinates
            # are not needed.
            if (Nz > 1):
                Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
                Zdim.Build(OutFile, Z)
                NoZ = False

        if (NoT):
            print("    Reading: {0:s} ({1:s})".format(InFname, Tname))
            T = InFile[Tname][...]
            Nt = len(T)
            Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
            Tdim.Build(OutFile, T)
            NoT = False

        # Read in input variable and form the domain average
        # Do one time step at a time in order to reduce memory usage
        for i in range(Nt):
            # Read in the var and determine if 2D (t,y,x) or 3D (t,z,y,x) field.
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
                    VarAvg[i,k] = CalcHorizAvg(np.squeeze(Var[i,:,:]), Select)

            if ((i % 10) == 0):
                print("    Processing time step: {0:d}".format(i))

        print("")

        InFile.close()

        # Write out average
        print("  Writing: {0:s} ({1:s})".format(OutFname, OutVname))
        Dset = h5u.DsetCoards(OutVname, Ndims, VarAvg.shape)
        if (Ndims == 2):
            Dset.Build(OutFile, VarAvg, Tdim)
        elif (Ndims == 3):
            Dset.Build(OutFile, VarAvg, Tdim, Zdim)

        print("")

    OutFile.close()
