#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigRce as conf
import Hdf5Utils as h5u
from scipy import signal as sig


Tstring = conf.SetTimeString()

SimList = [
#    'RCE_1km',
#    'RCE_1km_SM',
#    'RCE_1km_DM',
#    'RCE_1km_DP',
#    'RCE_1km_IPL2',
#    'RCE_1km_SM_IPL2',
#    'RCE_1km_DM_IPL2',
#    'RCE_1km_DP_IPL2',
    'RCE_3km_1mom',
    'RCE_3km_1mom_db',
    'RCE_3km_1mom_db_udef',
    'RCE_3km_1mom_db_rlongup',
    'RCE_3km_1mom_dm',
#    'RCE_3km_2mom',
    'RCE_3km_2mom_db',
    'RCE_3km_2mom_db_udef',
    'RCE_3km_2mom_db_rlongup',
    'RCE_3km_2mom_dm',
    'RCE_3km_2mom_dm_lrz',
    ]
Nsims = len(SimList)

VarList = [
    [ "dom_avg_pcprr",      "/pcprr"      ],
    [ "dom_avg_vint_vapor", "/vint_vapor" ],
    [ "dom_avg_top_lwup",   "/top_lwup"   ],
    ] 
Nvars = len(VarList)

# MAIN
for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating power spectra for simulation: {0:s}".format(Sim))
    print("")

    # Open the output file
    OutFname = "DIAGS/power_spectra_{0:s}.h5".format(Sim)
    OutFile = h5py.File(OutFname, mode='w')

    # Read input variables, calculate horizontal averages and write out results
    for ivar in range(Nvars):
        Fprefix = VarList[ivar][0]
        Vname   = VarList[ivar][1]

        # Open the input file

        InFname = "DIAGS/{0:s}_{1:s}.h5".format(Fprefix, Sim)
        InFile  = h5py.File(InFname, mode='r')

        # Read in input variable and form the temporal average using the
        # Tmin and Tmax values to select a range of time steps. Vars will
        # be of the form (t,z).
        print("  Reading: {0:s} ({1:s})".format(InFname, Vname))
        if (Vname == "/pcprr"):
            VAR = InFile[Vname][...] * 24.0 # convert mm/h to mm/day
        else:
            VAR = InFile[Vname][...]
        print("")

        # Calculate the power spectrum. Assume that the points in VAR are evenly spaced.
        VarFreq, VarPsd = sig.welch(VAR)
        Nf = len(VarFreq)

        # Write out average
        OutVname = "{0:s}_psd".format(Vname)
        print("  Writing: {0:s} ({1:s})".format(OutFname, OutVname))
        OutDset = h5u.BaseDset(OutVname, 1, [ Nf ])
        Dset = OutDset.Create(OutFile)
        Dset[...] = VarPsd

        OutVname = "{0:s}_freq".format(Vname)
        print("  Writing: {0:s} ({1:s})".format(OutFname, OutVname))
        OutDset = h5u.BaseDset(OutVname, 1, [ Nf ])
        Dset = OutDset.Create(OutFile)
        Dset[...] = VarFreq

        print("")

    InFile.close()
    OutFile.close()
