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
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_theta', '/ps_theta_hgrad', 'z' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_theta',  '/s_theta_hgrad',  'z' ],

    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_theta_e', '/ps_theta_e_hgrad', 'z' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_theta_e',  '/s_theta_e_hgrad',  'z' ],

    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_tempc',   '/ps_tempc_hgrad',   'z' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_tempc',    '/s_tempc_hgrad',    'z' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_tempc_p', '/ps_tempc_p_hgrad', 'p' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_tempc_p',  '/s_tempc_p_hgrad',  'p' ]
    ]
Nsets = len(PtrackList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Pname = '/p_coords'
Tname = '/t_coords'

# Need to do temp averaging between the near surface (110 m, 1000 mb)
# to the jet height (3200 m, 700 mb)
Zbot = 110  # m
Ztop = 3200

Pbot = 1000 # mb
Ptop = 700

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Creating PTRACK temperature gradients for simulation: {0:s}".format(Sim))
    print('')

    NoX = True
    NoY = True
    NoZ = True
    NoP = True
    NoT = True
    
    OutFname = "DIAGS/ptrack_tgrads_{0:s}.h5".format(Sim)
    Ofile = h5py.File(OutFname, mode='w')

    for iset in range(Nsets):
        InFname    = PtrackList[iset][0].replace("<SIM>", Sim)
        InVname    = PtrackList[iset][1]
        OutVname   = PtrackList[iset][2]
        VcoordType = PtrackList[iset][3]

        # Read input data. If this is the first set, read in the coordinates
        # and build the dimensions.
        print("  Reading {0:s} ({1:s})".format(InFname, InVname))
        Ifile = h5py.File(InFname, mode='r')
        Var = Ifile[InVname][...] # Var is (z,x), with the size of y-dim == 1

        if (NoX):
            print("    Reading {0:s} ({1:s})".format(InFname, Xname))
            X = Ifile[Xname][...]
            Nx = len(X)
            Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
            Xdim.Build(Ofile, X)

            # Create 1D delta X, repeat the final delta-x value to get
            # the length to match X.
            DeltaX = X[1:] - X[0:-1]
            DeltaX = np.append(DeltaX, DeltaX[-1])

            NoX = False

        if (NoY):
            print("    Reading {0:s} ({1:s})".format(InFname, Yname))
            Y = Ifile[Yname][...]
            Ny = len(Y)
            Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
            Ydim.Build(Ofile, Y)
            NoY = False

        if (VcoordType == 'z'):
            if (NoZ):
                print("    Reading {0:s} ({1:s})".format(InFname, Zname))
                Z = Ifile[Zname][...]
                Nz = len(Z)
                Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
                Zdim.Build(Ofile, Z)

                # Select Z range
                Select = np.where((Z >= Zbot) * (Z <= Ztop))
                Z1 = Select[0][0]
                Z2 = Select[0][-1]

                NoZ = False
        elif (VcoordType == 'p'):
            if (NoP):
                print("    Reading {0:s} ({1:s})".format(InFname, Pname))
                P = Ifile[Pname][...]
                Np = len(P)
                Pdim = h5u.DimCoards(Pname, 1, [ Np ], 'p')
                Pdim.Build(Ofile, P)

                # Select P range
                Select = np.where((P <= Pbot) * (P >= Ptop))
                P1 = Select[0][0]
                P2 = Select[0][-1]

                NoP = False

        if (NoT):
            print("    Reading {0:s} ({1:s})".format(InFname, Tname))
            T = Ifile[Tname][...]
            Nt = len(T)
            Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
            Tdim.Build(Ofile, T)
            NoT = False

        Ifile.close()

        # Do vertical average of layer between bottom and top layer
        if (VcoordType == 'z'):
            Tbot = np.squeeze(Var[Z1,:])
            Ttop = np.squeeze(Var[Z2,:])
        elif (VcoordType == 'p'):
            Tbot = np.squeeze(Var[P1,:])
            Ttop = np.squeeze(Var[P2,:])
        Tbar = (Tbot + Ttop) * 0.5

        # Calc gradient: DeltaT/DeltaX
        DeltaT = Tbar[1:] - Tbar[:-1]
        DeltaT = np.append(DeltaT, DeltaT[-1])

        Var = DeltaT / DeltaX
        VarSmooth = pltu.SmoothLine(Var, 5)
         
        # Write out gradient, plus smoothed version of gradient
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        VarDset = h5u.DsetCoards(OutVname, 1, [ Nx ])
        VarDset.Build(Ofile, Var, Xdim)

        Vname = "{0:s}_smooth".format(OutVname)
        print("  Writing {0:s} ({1:s})".format(OutFname, Vname))
        VarDset = h5u.DsetCoards(Vname, 1, [ Nx ])
        VarDset.Build(Ofile, VarSmooth, Xdim)
        print('')

    Ofile.close()
