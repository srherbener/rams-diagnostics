#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'  # vertical coords: height
Pname = '/p_coords'  # vertical coords: pressure
Tname = '/t_coords'


for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Removing vortex from horizontal winds: {0:s}".format(Sim))
    print('')

    InFname = "HDF5/{0:s}/HDF5/u_lite-{0:s}-AS-2006-08-20-120000-g3.h5".format(Sim, Sim)
    Uvname = '/u'
    Vvname = '/v'




#    NoX = True
#    NoY = True
#    NoZ = True
#    NoP = True
#    NoT = True
#
#    OutFname = "DIAGS/ptrack_avgs_{0:s}.h5".format(Sim)
#    Ofile = h5py.File(OutFname, mode='w')
#
#    for iset in range(Nsets):
#        InFname    = PtrackList[iset][0].replace("<SIM>", Sim)
#        InVname    = PtrackList[iset][1]
#        AvgPeriod  = PtrackList[iset][2]
#        OutVname   = PtrackList[iset][3]
#        VcoordType = PtrackList[iset][4]
#    
#        # Read input data. If this is the first set, read in the coordinates
#        # and build the dimensions.
#        print("  Reading {0:s} ({1:s})".format(InFname, InVname))
#        print("    Averaging time period: {0:s}".format(AvgPeriod))
#        print("    Vertical coordinate type: {0:s}".format(VcoordType))
#        Ifile = h5py.File(InFname, mode='r')
#        Var = np.squeeze(Ifile[InVname][...]) # Y-dim is size 1, reduce Var to (t,z,x)
#
#        if (NoX):
#            print("    Reading {0:s} ({1:s})".format(InFname, Xname))
#            X = Ifile[Xname][...]
#            Nx = len(X)
#            Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
#            Xdim.Create(Ofile, X)
#            NoX = False
#
#        if (NoY):
#            print("    Reading {0:s} ({1:s})".format(InFname, Yname))
#            Y = Ifile[Yname][...]
#            Ny = len(Y)
#            Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
#            Ydim.Create(Ofile, Y)
#            NoY = False
#
#        if (VcoordType == 'z'):
#            if (NoZ):
#                print("    Reading {0:s} ({1:s})".format(InFname, Zname))
#                Z = Ifile[Zname][...]
#                Nz = len(Z)
#                Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
#                Zdim.Create(Ofile, Z)
#                NoZ = False
#        elif (VcoordType == 'p'):
#            if (NoP):
#                print("    Reading {0:s} ({1:s})".format(InFname, Zname))
#                P = Ifile[Zname][...]
#                Np = len(P)
#                Pdim = h5u.DimCoards(Pname, 1, [ Np ], 'p')
#                Pdim.Create(Ofile, P)
#                NoP = False
#
#        if (NoT):
#            print("    Reading {0:s} ({1:s})".format(InFname, Tname))
#            T = Ifile[Tname][...]
#            SimT = T / 3600.0 - 42
#            Nt = len(T)
#            Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
#            Tdim.Create(Ofile, T)
#
#            # Find indices of PSAP and SAP time periods
#            SimT = T / 3600.0 - 42
#
#            Select = np.where((SimT >= PsapStart) * (SimT <= PsapEnd))
#            PS_T1 = Select[0][0]
#            PS_T2 = Select[0][-1]
#
#            Select = np.where((SimT >= SapStart) * (SimT <= SapEnd))
#            S_T1 = Select[0][0]
#            S_T2 = Select[0][-1]
#
#            NoT = False
#
#        Ifile.close()
#
#        # Want to reduce time dimensions by taking mean value across
#        # the time period.
#        if (AvgPeriod == 'pre_sal'):
#            AvgVar = np.squeeze(np.mean(Var[PS_T1:PS_T2,:,:], axis=0))
#        elif (AvgPeriod == 'sal'):
#            AvgVar = np.squeeze(np.mean(Var[S_T1:S_T2,:,:], axis=0))
#
#        # Write out PSAP and SAP averged vars
#        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
#        if (VcoordType == 'z'):
#            VarDset = h5u.DsetCoards(OutVname, 2, [ Nz, Nx ])
#            VarDset.Create(Ofile, AvgVar, Zdim, Xdim)
#        elif (VcoordType == 'p'):
#            VarDset = h5u.DsetCoards(OutVname, 2, [ Np, Nx ])
#            VarDset.Create(Ofile, AvgVar, Pdim, Xdim)
#        print('')
#
#    Ofile.close()
