#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))

import h5py
import numpy as np
from matplotlib import pyplot as plt
import ConfigRce as conf
import Hdf5Utils as h5u
import PlotUtils as plu


Tstring = conf.SetTimeString()

SimList = [
#    'RCE_3km_1mom',
#    'RCE_3km_1mom_db',
    'RCE_3km_1mom_db_udef',
    'RCE_3km_1mom_db_rlongup',
#    'RCE_3km_1mom_dm',
#    'RCE_3km_2mom',
#    'RCE_3km_2mom_db',
    'RCE_3km_2mom_db_udef',
    'RCE_3km_2mom_db_rlongup',
#    'RCE_3km_2mom_dm',
#    'RCE_3km_2mom_dm_lrz',
    ]
Nsims = len(SimList)

VarList = [
    #  <in_file_prefix> <in_var_name> <increment_between_frames> <Cmin> <Cmax>
    #
    [ "pcprr",      "/pcprr",        4, "Prate", 1e-3, 1e2, 11, "nipy_spectral", "log" ],
    [ "vint_tcond", "/vertint_orig", 4, "Tcond", 1e-3, 1e2, 11, "nipy_spectral", "log" ],

    ] 
Nvars = len(VarList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

InFileTemplate = "SIMS/<SIM>/HDF5/<FPREFIX>-<SIM>-LS-2012-01-01-000000-g1.h5"
OutDirTemplate = "Animations/<SIM>/<OUTVAR>"

Cfilled = True

for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating animation frames for simulation: {0:s}".format(Sim))
    print("")

    for ivar in range(Nvars):
        InFprefix = VarList[ivar][0]
        InVname   = VarList[ivar][1]
        FrameInc  = VarList[ivar][2]
        OutVname  = VarList[ivar][3]
        Cmin      = VarList[ivar][4]
        Cmax      = VarList[ivar][5]
        Cnum      = VarList[ivar][6]
        Cmap      = VarList[ivar][7]
        Ctype     = VarList[ivar][8]

        InFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", InFprefix)
        InFile  = h5py.File(InFname, mode='r')
        print("  Reading: {0:s} ({1:s})".format(InFname, InVname))
        print("    Reading: {0:s} ({1:s})".format(InFname, Xname))
        print("    Reading: {0:s} ({1:s})".format(InFname, Yname))
        print("    Reading: {0:s} ({1:s})".format(InFname, Tname))

        X = InFile[Xname][...]
        Y = InFile[Yname][...]
        T = InFile[Tname][...]

        Nt = len(T)
        OutNt = (Nt // FrameInc) + 1 # // is integer divide

        OutDir = OutDirTemplate.replace("<SIM>", Sim).replace("<OUTVAR>", OutVname)
        if (not os.path.isdir(OutDir)):
            os.makedirs(OutDir)
        print("  Writing {0:d} frames into: {1:s}".format(OutNt, OutDir))

        Count = 0
        for it in range(0,Nt,FrameInc):
            Count = Count + 1
            SimTime = T[it] / 86400.0 # convert to days

            # read in the variable, it will keep track of the frame increment
            VAR = InFile[InVname][it,...]

            OutFname = "{0:s}/{1:s}_{2:s}_{3:04d}.png".format(OutDir, OutVname, Sim, Count)
            print("    {0:s}".format(OutFname))

            Fig = plt.figure()
            Xaxis = plu.AxisConfig('x', [ X[0], X[-1] ], "Longitude")
            Xaxis.fontsize = 14
            Yaxis = plu.AxisConfig('y', [ Y[0], Y[-1] ], "Latitude")
            Yaxis.fontsize = 14

            Cspecs = plu.ContourConfig(Cmin, Cmax, Cnum, Cmap, Cfilled, Ctype)

            TitleString = "{0:s}: {1:s}, SimTime: {2:7.2f} (d)".format(Sim, OutVname, SimTime)
            Ptitle = plu.TitleConfig("", TitleString)
            Ptitle.fontsize = 20

            plu.PlotContour(Fig.gca(), X, Y, VAR, Ptitle, Xaxis, Yaxis, Cspecs)

            Fig.savefig(OutFname)
            plt.close()

        InFile.close()
        print("")

