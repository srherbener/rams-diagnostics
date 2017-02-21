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
RhoAir = conf.SetRhoAir()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

PtrackList = [
    [ 'XsectionData/ptrack_theta_<SIM>.h5', '/theta',  'none',  'pre_sal', '/ps_theta' ],
    [ 'XsectionData/ptrack_theta_<SIM>.h5', '/theta',  'none',  'sal',     '/s_theta'  ],
    [ 'XsectionData/ptrack_theta_<SIM>.h5', '/theta',  'hgrad', 'pre_sal', '/ps_theta_hgrad' ],
    [ 'XsectionData/ptrack_theta_<SIM>.h5', '/theta',  'hgrad', 'sal',     '/s_theta_hgrad'  ],

    [ 'XsectionData/ptrack_theta_e_<SIM>.h5', '/theta_e',  'none',  'pre_sal', '/ps_theta_e' ],
    [ 'XsectionData/ptrack_theta_e_<SIM>.h5', '/theta_e',  'none',  'sal',     '/s_theta_e'  ],
    [ 'XsectionData/ptrack_theta_e_<SIM>.h5', '/theta_e',  'hgrad', 'pre_sal', '/ps_theta_e_hgrad' ],
    [ 'XsectionData/ptrack_theta_e_<SIM>.h5', '/theta_e',  'hgrad', 'sal',     '/s_theta_e_hgrad'  ],

    [ 'XsectionData/ptrack_tempc_<SIM>.h5', '/tempc',  'none',  'pre_sal', '/ps_tempc' ],
    [ 'XsectionData/ptrack_tempc_<SIM>.h5', '/tempc',  'none',  'sal',     '/s_tempc'  ],
    [ 'XsectionData/ptrack_tempc_<SIM>.h5', '/tempc',  'hgrad', 'pre_sal', '/ps_tempc_hgrad' ],
    [ 'XsectionData/ptrack_tempc_<SIM>.h5', '/tempc',  'hgrad', 'sal',     '/s_tempc_hgrad'  ]
    ]
Nsets = len(PtrackList)

PsapStart = 10  # PSAP, sim time in hours
PsapEnd   = 30

SapStart = 40  # SAP, sim time in hours
SapEnd   = 60

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Creating PTRACK averages for simulation: {0:s}".format(Sim))
    print('')

    OutFname = "DIAGS/ptrack_avgs_{0:s}.h5".format(Sim)
    Ofile = h5py.File(OutFname, mode='w')

    for iset in range(Nsets):
        InFname   = PtrackList[iset][0].replace("<SIM>", Sim)
        InVname   = PtrackList[iset][1]
        InMod     = PtrackList[iset][2]
        AvgPeriod = PtrackList[iset][3]
        OutVname  = PtrackList[iset][4]
    
        # Read input data. If this is the first set, read in the coordinates
        # and build the dimensions.
        print("  Reading {0:s} ({1:s})".format(InFname, InVname))
        print("    Input data modification: {0:s}".format(InMod))
        print("    Averaging time period: {0:s}".format(AvgPeriod))
        Ifile = h5py.File(InFname, mode='r')
        Var = Ifile[InVname][...] # Var is (t,z,y,x), with the size of y-dim == 1
        Var = np.squeeze(Var)     # Var is now (t,z,x)

        if (iset == 0):
            X = Ifile[Xname][...]
            Y = Ifile[Yname][...]
            Z = Ifile[Zname][...]
            T = Ifile[Tname][...]
            SimT = T / 3600.0 - 42

            Nx = len(X)
            Ny = len(Y)
            Nz = len(Z)
            Nt = len(T)

            Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
            Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
            Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
            Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)

            Xdim.Create(Ofile, X)
            Ydim.Create(Ofile, Y)
            Zdim.Create(Ofile, Z)
            Tdim.Create(Ofile, T)

            # Find indices of PSAP and SAP time periods
            SimT = T / 3600.0 - 42

            Select = np.where((SimT >= PsapStart) * (SimT <= PsapEnd))
            PS_T1 = Select[0][0]
            PS_T2 = Select[0][-1]

            Select = np.where((SimT >= SapStart) * (SimT <= SapEnd))
            S_T1 = Select[0][0]
            S_T2 = Select[0][-1]

            # Create 2D gridded X and Y values for doing spatial gradients
            # Repeat the 2D grids along the time dimension
            # Create delta X by taking differences of adjacent cells
            # along the x-dimension. Repeat the final x value to get shape
            # back to [Nt,Nz,Nx].
            X2D, Z2D = np.meshgrid(X, Z, indexing='xy')
            X3D = np.zeros([Nt, Nz, Nx])
            for i in range(Nt):
                X3D[i,:,:] = X2D
            DeltaX = X3D[:,:,1:] - X3D[:,:,0:-1]
            DeltaX = np.append(DeltaX, DeltaX[:,:,[-1]], axis=2)

        Ifile.close()

        # If requested, set Var to horizontal gradient of Var
        if (InMod == 'hgrad'):
            DeltaVar = Var[:,:,1:] - Var[:,:,0:-1]
            DeltaVar = np.append(DeltaVar, DeltaVar[:,:,[-1]], axis=2)
            Var = DeltaVar / DeltaX
          
        # Want to reduce time dimensions by taking mean value across
        # the time period.
        if (AvgPeriod == 'pre_sal'):
            AvgVar = np.squeeze(np.mean(Var[PS_T1:PS_T2,:,:], axis=0))
        elif (AvgPeriod == 'sal'):
            AvgVar = np.squeeze(np.mean(Var[S_T1:S_T2,:,:], axis=0))

        # Write out PSAP and SAP averged vars
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        VarDset = h5u.DsetCoards(OutVname, 2, [ Nz, Nx ])
        VarDset.Create(Ofile, AvgVar, Zdim, Xdim)
        print('')

    Ofile.close()
