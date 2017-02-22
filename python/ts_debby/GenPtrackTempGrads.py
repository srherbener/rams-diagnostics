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

PtrackList = [
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_theta', '/ps_theta_hgrad' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_theta',  '/s_theta_hgrad'  ],

    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_theta_e', '/ps_theta_e_hgrad' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_theta_e',  '/s_theta_e_hgrad'  ],

    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_tempc', '/ps_tempc_hgrad' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_tempc',  '/s_tempc_hgrad'  ]
    ]
Nsets = len(PtrackList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

Zmin = 110  # height near surface (1000 mb)
Zmax = 3500 # height of mid-level jet (658 mb)

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Creating PTRACK temperature gradients for simulation: {0:s}".format(Sim))
    print('')

    OutFname = "DIAGS/ptrack_tgrads_{0:s}.h5".format(Sim)
    Ofile = h5py.File(OutFname, mode='w')

    for iset in range(Nsets):
        InFname   = PtrackList[iset][0].replace("<SIM>", Sim)
        InVname   = PtrackList[iset][1]
        OutVname  = PtrackList[iset][2]
    
        # Read input data. If this is the first set, read in the coordinates
        # and build the dimensions.
        print("  Reading {0:s} ({1:s})".format(InFname, InVname))
        Ifile = h5py.File(InFname, mode='r')
        Var = Ifile[InVname][...] # Var is (z,x), with the size of y-dim == 1

        if (iset == 0):
            X = Ifile[Xname][...]
            Y = Ifile[Yname][...]
            Z = Ifile[Zname][...]
            T = Ifile[Tname][...]

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

            # Create 1D delta X, repeat the final delta-x value to get
            # the length to match X.
            DeltaX = X[1:] - X[0:-1]
            DeltaX = np.append(DeltaX, DeltaX[-1])

            # Select Z range
            Select = np.where((Z >= Zmin) * (Z <= Zmax))
            Z1 = Select[0][0]
            Z2 = Select[0][-1]

        Ifile.close()

        # Do vertical average of layer between Zmin and Zmax
        Tbot = np.squeeze(Var[Z1,:])
        Ttop = np.squeeze(Var[Z2,:])
        Tbar = (Tbot + Ttop) * 0.5

        # Calc gradient: DeltaT/DeltaX
        DeltaT = Tbar[1:] - Tbar[:-1]
        DeltaT = np.append(DeltaT, DeltaT[-1])

        Var = DeltaT / DeltaX
        VarSmooth = pltu.SmoothLine(Var, 5)
         
        # Write out gradient, plus smoothed version of gradient
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        VarDset = h5u.DsetCoards(OutVname, 1, [ Nx ])
        VarDset.Create(Ofile, Var, Xdim)

        Vname = "{0:s}_smooth".format(OutVname)
        print("  Writing {0:s} ({1:s})".format(OutFname, Vname))
        VarDset = h5u.DsetCoards(Vname, 1, [ Nx ])
        VarDset.Create(Ofile, VarSmooth, Xdim)
        print('')

    Ofile.close()
