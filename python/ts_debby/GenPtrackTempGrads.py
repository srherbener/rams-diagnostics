#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import scipy.integrate as integ
import pandas as pd
import ConfigTsd as conf
import Hdf5Utils as h5u

#######################################################################################
# Functions
#######################################################################################

#####################################
# SmoothLine
#
# Remove noise from a linear data series.
#
def SmoothLine(Var):
    # Use the ewma (exponentially weighted moving average) from the
    # scipy package to smooth a data series. Apply both forward and
    # backward in order to eliminate phase delay.

    FiltCom = 10 # Center of mass spec for ewma function. alpha = 1 / ( 1 + FiltCom)

    Forward  = pd.ewma(Var, com=FiltCom)
    Backward = pd.ewma(Var[::-1], com=FiltCom)
    VarSmooth = np.vstack((Forward, Backward[::-1]))
    VarSmooth = np.squeeze(np.mean(VarSmooth, axis=0))

    return VarSmooth

#######################################################################################
# Main
#######################################################################################

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
#    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_tempc',   'z' ],
#    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_tempc',    'z' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/ps_tempc_nv_p', 'p' ],
    [ 'DIAGS/ptrack_avgs_<SIM>.h5', '/s_tempc_nv_p',  'p' ],
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

FiltLen = 5

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
        VcoordType = PtrackList[iset][2]

        # Read input data. If this is the first set, read in the coordinates
        # and build the dimensions.
        print("  Reading {0:s} ({1:s})".format(InFname, InVname))
        Ifile = h5py.File(InFname, mode='r')
        InVar = Ifile[InVname][...] # InVar is (z,x), with the size of y-dim == 1

        if (NoX):
            print("    Reading {0:s} ({1:s})".format(InFname, Xname))
            X = Ifile[Xname][...]
            Nx = len(X)
            Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
            Xdim.Build(Ofile, X)

            # Create delta X: dx[i] = x[i+1] - x[i]
            # then repeat the last dx[i] to make delta X length
            # match that of X
            DeltaX = np.diff(X)
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

        # Record items that depend on VcoordType
        if (VcoordType == 'z'):
            Nvert = Nz

            Height = Z[Z2] - Z[Z1]
            Hselect = Z[Z1:Z2]

            Tselect = np.squeeze(InVar[Z1:Z2,:])

            VertDim = Zdim
        elif (VcoordType == 'p'):
            Nvert = Np

            Height = P[P2] - P[P1]
            Hselect = P[P1:P2]

            Tselect = np.squeeze(InVar[P1:P2,:])

            VertDim = Pdim

        # Smooth the input temperature field in the horizontal.
        InVarSmooth = np.zeros((Nvert,Nx))
        for ivert in range(Nvert):
            InVarSmooth[ivert,:] = SmoothLine(np.squeeze(InVar[ivert,:]))

        # Take the horizontal gradient of the 2D field.
        #   Create a 2D delta X array
        DeltaX_2D = np.tile(DeltaX.reshape((1,Nx)), (Nvert,1))
        InVarGrad = np.gradient(InVarSmooth, DeltaX_2D, axis=1)
        InVarGradSmooth = np.zeros((Nvert,Nx))
        for ivert in range(Nvert):
            InVarGradSmooth[ivert,:] = SmoothLine(np.squeeze(InVarGrad[ivert,:]))

        # Do vertical average of layer between bottom and top layer
        # Use Simpson's Rule integration to the do the averaging
        Hselect = Hselect.reshape((len(Hselect),1))
        Hselect = np.tile(Hselect, (1, Nx))
        Tbar = np.squeeze(integ.simps(Tselect, Hselect, axis=0)) / Height

        # Smooth Tbar in order to mitigate noise in the gradient
        Tbar = SmoothLine(Tbar)

        # Gradient of Tbar (mean temperature of layer)
        Tgrad = np.gradient(Tbar, DeltaX)
        TgradSmooth = SmoothLine(Tgrad)
         
        # Output:
        #  original field
        #  smoothed field
        #  gradient of field
        #  smoothed gradient of field
        #  Tbar
        #  gradient based on Tbar
        #  smoothed gradient based on Tbar
        print("  Writing {0:s} ({1:s})".format(OutFname, InVname))
        TempDset = h5u.DsetCoards(InVname, 2, [ Nvert, Nx ])
        TempDset.Build(Ofile, InVar, VertDim, Xdim)

        OutVname = "{0:s}_smooth".format(InVname)
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        TempDset = h5u.DsetCoards(OutVname, 2, [ Nvert, Nx ])
        TempDset.Build(Ofile, InVarSmooth, VertDim, Xdim)

        OutVname = "{0:s}_hgrad".format(InVname)
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        TempDset = h5u.DsetCoards(OutVname, 2, [ Nvert, Nx ])
        TempDset.Build(Ofile, InVarGrad, VertDim, Xdim)

        OutVname = "{0:s}_hgrad_smooth".format(InVname)
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        TempDset = h5u.DsetCoards(OutVname, 2, [ Nvert, Nx ])
        TempDset.Build(Ofile, InVarGradSmooth, VertDim, Xdim)

        OutVname = "{0:s}_bar".format(InVname)
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        TbarDset = h5u.DsetCoards(OutVname, 1, [ Nx ])
        TbarDset.Build(Ofile, Tbar, Xdim)

        OutVname = "{0:s}_bar_hgrad".format(InVname)
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        TgradDset = h5u.DsetCoards(OutVname, 1, [ Nx ])
        TgradDset.Build(Ofile, Tgrad, Xdim)

        OutVname = "{0:s}_bar_hgrad_smooth".format(InVname)
        print("  Writing {0:s} ({1:s})".format(OutFname, OutVname))
        TgradSmoothDset = h5u.DsetCoards(OutVname, 1, [ Nx ])
        TgradSmoothDset.Build(Ofile, TgradSmooth, Xdim)
        print('')

    Ofile.close()
