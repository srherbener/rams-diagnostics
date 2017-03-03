#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u
import PlotUtils as pltu


Tstring = conf.SetTimeString()
RhoAir = conf.SetRhoAir()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

Xname  = '/x_coords'
XLname = '/xl_coords'
Yname  = '/y_coords'
Zname  = '/z_coords'
Pname  = '/p_coords'
Tname  = '/t_coords'

PreSalStart = 10 # sim time in hours
PreSalEnd   = 30
SalStart    = 40
SalEnd      = 60

Pdx = 1.000
Pdy = 2.593
Pslope = Pdy / Pdx
Pangle = np.arctan2(Pdy, Pdx)

PtrackList = [
    [ 'XsectionData/ptrack_u_<SIM>.h5',   '/u', 'XsectionData/ptrack_v_<SIM>.h5',   '/v', 'u',   'v',    'x', 'z' ],
    [ 'XsectionData/ptrack_u_p_<SIM>.h5', '/u', 'XsectionData/ptrack_v_p_<SIM>.h5', '/v', 'u_p', 'v_p',  'x', 'p' ],

    [ 'XsectionData/ptrack_u_nv_lite_<SIM>.h5',   '/u', 'XsectionData/ptrack_v_nv_lite_<SIM>.h5',   '/v', 'u_nv_lite',   'v_nv_lite',    'xl', 'z' ],
    [ 'XsectionData/ptrack_u_nv_lite_p_<SIM>.h5', '/u', 'XsectionData/ptrack_v_nv_lite_p_<SIM>.h5', '/v', 'u_nv_lite_p', 'v_nv_lite_p',  'xl', 'p' ]
    ]
Nsets = len(PtrackList)

print("Translating horizontal velocity vectors to ptrack axes")
print("  Ptrack dx, dy: {0:.3f}, {1:.3f}".format(Pdx, Pdy))
print("  Ptrack slope: {0:.3f}".format(Pslope))
print("  Ptrack angle: {0:.3f}".format(Pangle))
print("")

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Simulation: {0:s}".format(Sim))
    print('')

    NoX  = True
    NoXL = True
    NoY  = True
    NoZ  = True
    NoP  = True
    NoT  = True
    
    OutFname = "DIAGS/ptrack_hvel_{0:s}.h5".format(Sim)
    Ofile = h5py.File(OutFname, mode='w')

    for iset in range(Nsets):
        InUfname   = PtrackList[iset][0].replace("<SIM>", Sim)
        InUvname   = PtrackList[iset][1]
        InVfname   = PtrackList[iset][2].replace("<SIM>", Sim)
        InVvname   = PtrackList[iset][3]
        OutUvname  = PtrackList[iset][4]
        OutVvname  = PtrackList[iset][5]
        XcoordType = PtrackList[iset][6]
        VcoordType = PtrackList[iset][7]

        # Read input data. If this is the first set, read in the coordinates
        # and build the dimensions.
        print("  Reading {0:s} ({1:s})".format(InUfname, InUvname))
        Ufile = h5py.File(InUfname, mode='r')
        U = np.squeeze(Ufile[InUvname][...]) # U is (t,z,x)

        print("  Reading {0:s} ({1:s})".format(InVfname, InVvname))
        Vfile = h5py.File(InVfname, mode='r')
        V = np.squeeze(Vfile[InVvname][...]) # V is (t,z,x)

        if (XcoordType == 'x'):
            if (NoX):
                print("    Reading {0:s} ({1:s})".format(InUfname, Xname))
                X = Ufile[Xname][...]
                Nx = len(X)
                Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
                Xdim.Build(Ofile, X)
                NoX = False
        elif (XcoordType == 'xl'):
            if (NoXL):
                print("    Reading {0:s} ({1:s})".format(InUfname, Xname))
                XL = Ufile[Xname][...]
                Nxl = len(XL)
                XLdim = h5u.DimCoards(XLname, 1, [ Nxl ], 'x')
                XLdim.Build(Ofile, XL)
                NoXL = False

        if (NoY):
            print("    Reading {0:s} ({1:s})".format(InUfname, Yname))
            Y = Ufile[Yname][...]
            Ny = len(Y)
            Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
            Ydim.Build(Ofile, Y)
            NoY = False

        if (VcoordType == 'z'):
            if (NoZ):
                print("    Reading {0:s} ({1:s})".format(InUfname, Zname))
                Z = Ufile[Zname][...]
                Nz = len(Z)
                Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
                Zdim.Build(Ofile, Z)
                NoZ = False
        elif (VcoordType == 'p'):
            if (NoP):
                print("    Reading {0:s} ({1:s})".format(InUfname, Zname))
                P = Ufile[Zname][...]
                Np = len(P)
                Pdim = h5u.DimCoards(Pname, 1, [ Np ], 'p')
                Pdim.Build(Ofile, P)
                NoP = False

        if (NoT):
            print("    Reading {0:s} ({1:s})".format(InUfname, Tname))
            T = Ufile[Tname][...]
            Nt = len(T)
            Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
            Tdim.Build(Ofile, T)
            NoT = False

            # Select time ranges for PSAP and SAP
            SIM_T = T / 3600.0 - 42.0 # convert seconds to hours

            Select = np.where((SIM_T >= PreSalStart) * (SIM_T <= PreSalEnd))
            PS_T1 = Select[0][0]
            PS_T2 = Select[0][-1]

            Select = np.where((SIM_T >= SalStart) * (SIM_T <= SalEnd))
            S_T1 = Select[0][0]
            S_T2 = Select[0][-1]

        Ufile.close()
        Vfile.close()

        # Convert to polar coordinates (relative to domain axes)
        Mag   = np.sqrt(np.square(U) + np.square(V))
        Angle = np.arctan2(V, U)

        # Translate to ptrack axes
        Angle = Angle - Pangle

        # Convert to Cartesian coordinates (relative to ptrack axes)
        Uptrack = Mag * np.cos(Angle)
        Vptrack = Mag * np.sin(Angle)

        # Do the temporal averaging
        PS_Uptrack = np.squeeze(np.mean(Uptrack[PS_T1:PS_T2,:,:], axis=0)) # result is (z,x)
        PS_Vptrack = np.squeeze(np.mean(Vptrack[PS_T1:PS_T2,:,:], axis=0))

        S_Uptrack = np.squeeze(np.mean(Uptrack[S_T1:S_T2,:,:], axis=0))
        S_Vptrack = np.squeeze(np.mean(Vptrack[S_T1:S_T2,:,:], axis=0))

        # Write out translated u, v, and temporally averaged u, v
        if (XcoordType == 'x'):
            HorizDim = Xdim
        elif (XcoordType == 'xl'):
            HorizDim = XLdim

        if (VcoordType == 'z'):
            VertDim = Zdim
        elif (VcoordType == 'p'):
            VertDim = Pdim

        Vname = "/{0:s}".format(OutUvname)
        print("  Writing {0:s} ({1:s})".format(OutFname, Vname))
        VarDset = h5u.DsetCoards(Vname, Uptrack.ndim, Uptrack.shape)
        VarDset.Build(Ofile, Uptrack, Tdim, VertDim, HorizDim)

        Vname = "/{0:s}".format(OutVvname)
        print("  Writing {0:s} ({1:s})".format(OutFname, Vname))
        VarDset = h5u.DsetCoards(Vname, Vptrack.ndim, Vptrack.shape)
        VarDset.Build(Ofile, Vptrack, Tdim, VertDim, HorizDim)


        Vname = "/ps_{0:s}".format(OutUvname)
        print("  Writing {0:s} ({1:s})".format(OutFname, Vname))
        VarDset = h5u.DsetCoards(Vname, PS_Uptrack.ndim, PS_Uptrack.shape)
        VarDset.Build(Ofile, PS_Uptrack, VertDim, HorizDim)

        Vname = "/ps_{0:s}".format(OutVvname)
        print("  Writing {0:s} ({1:s})".format(OutFname, Vname))
        VarDset = h5u.DsetCoards(Vname, PS_Vptrack.ndim, PS_Vptrack.shape)
        VarDset.Build(Ofile, PS_Vptrack, VertDim, HorizDim)

        Vname = "/s_{0:s}".format(OutUvname)
        print("  Writing {0:s} ({1:s})".format(OutFname, Vname))
        VarDset = h5u.DsetCoards(Vname, S_Uptrack.ndim, S_Uptrack.shape)
        VarDset.Build(Ofile, S_Uptrack, VertDim, HorizDim)

        Vname = "/s_{0:s}".format(OutVvname)
        print("  Writing {0:s} ({1:s})".format(OutFname, Vname))
        VarDset = h5u.DsetCoards(Vname, S_Vptrack.ndim, S_Vptrack.shape)
        VarDset.Build(Ofile, S_Vptrack, VertDim, HorizDim)

        print('')

    Ofile.close()
