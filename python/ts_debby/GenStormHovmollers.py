#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
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

HovmollerList = [
    [ 'DIAGS/hist_meas_az_vapor_<SIM>.h5', '/all_core_vapor_ts', '/all_core_vapor' ],
    [ 'DIAGS/hist_meas_az_vapor_<SIM>.h5', '/all_rb_vapor_ts',   '/all_rb_vapor'   ],

    [ 'DIAGS/hist_meas_az_theta_<SIM>.h5', '/all_core_theta_ts', '/all_core_theta' ],
    [ 'DIAGS/hist_meas_az_theta_<SIM>.h5', '/all_rb_theta_ts',   '/all_rb_theta'   ],

    [ 'DIAGS/hist_meas_az_w_<SIM>.h5', '/all_core_updraft_ts', '/all_core_updraft' ],
    [ 'DIAGS/hist_meas_az_w_<SIM>.h5', '/all_rb_updraft_ts',   '/all_rb_updraft'   ],
    [ 'DIAGS/hist_meas_az_w_<SIM>.h5', '/all_core_dndraft_ts', '/all_core_dndraft' ],
    [ 'DIAGS/hist_meas_az_w_<SIM>.h5', '/all_rb_dndraft_ts',   '/all_rb_dndraft'   ]

    ]
Nsets = len(HovmollerList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Creating Hovmoller data for simulation: {0:s}".format(Sim))
    print('')

    OutFname = "DIAGS/storm_hovs_{0:s}.h5".format(Sim)
    Ofile = h5py.File(OutFname, mode='w')

    for iset in range(Nsets):
        InFname  = HovmollerList[iset][0].replace("<SIM>", Sim)
        InVname  = HovmollerList[iset][1]
        OutVname = HovmollerList[iset][2]
    
        # Read input data. If this is the first set, read in the coordinates
        # and build the dimensions.
        print("  Reading {0:s} ({1:s})".format(InFname, InVname))
        Ifile = h5py.File(InFname, mode='r')
        Var = Ifile[InVname][...]
        Var = Var.transpose() # transpose since 2D contour plot routines expect
                              # data to be in the form (y,x), or for these
                              # hovmollers (z,t)
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

        Ifile.close()

        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        VarDset = h5u.DsetCoards(OutVname, 2, [ Nz, Nt ])
        VarDset.Create(Ofile, Var, Zdim, Tdim)
        print('')

    Ofile.close()
