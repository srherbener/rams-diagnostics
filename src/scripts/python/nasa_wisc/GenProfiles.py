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
    'RCE_3km_2mom',
    'RCE_3km_2mom_db',
    'RCE_3km_2mom_db_udef',
    'RCE_3km_2mom_db_rlongup',
    'RCE_3km_2mom_dm',
    'RCE_3km_2mom_dm_lrz',
    ]
Nsims = len(SimList)

# Start and end times are time step indexes. Each time step is one hour and
# we want the final 20 days.
#   T1 = -481   480 timesteps from the end (20 X 24 hours)
#   T2 = -1     end
VarList = [
    [ "dom_avg_theta",      "/theta",      -481, -1, 1.0, 0.0 ],
    [ "dom_avg_total_cond", "/total_cond", -481, -1, 1.0, 0.0 ],
    [ "dom_avg_vapor",      "/vapor",      -481, -1, 1.0, 0.0 ],
    [ "dom_avg_tempk",      "/tempk",      -481, -1, 1.0, 0.0 ],
    [ "dom_avg_relhum",     "/relhum",     -481, -1, 1.0, 0.0 ],
    ] 
Nvars = len(VarList)

Zname = '/z_coords'

# MAIN
for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating profiles for simulation: {0:s}".format(Sim))
    print("")

    # Open the output file
    OutFname = "DIAGS/profiles_{0:s}.h5".format(Sim)
    OutFile = h5py.File(OutFname, mode='w')

    # Read input variables, calculate horizontal averages and write out results
    for ivar in range(Nvars):
        Fprefix = VarList[ivar][0]
        Vname   = VarList[ivar][1]
        T1      = VarList[ivar][2]
        T2      = VarList[ivar][3]
        Scale   = VarList[ivar][4]
        Offset  = VarList[ivar][5]

        # Open the input file

        InFname = "DIAGS/{0:s}_{1:s}.h5".format(Fprefix, Sim)
        InFile  = h5py.File(InFname, mode='r')

        # Read in Z coordinates and build the Z dimension in the output file
        if (ivar == 0):
            Z = InFile[Zname][...]
            Nz = len(Z)
            Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
            Zdim.Build(OutFile, Z)

        # Read in input variable and form the temporal average using the
        # Tmin and Tmax values to select a range of time steps. Vars will
        # be of the form (t,z).
        print("  Reading: {0:s} ({1:s}, {2:e}, {3:e})".format(InFname, Vname, Scale, Offset))
        Var = (np.mean(InFile[Vname][T1:T2,...], axis=0) * Scale) + Offset

        # Write out average
        OutVname = "{0:s}_prof".format(Vname)
        print("  Writing: {0:s} ({1:s})".format(OutFname, OutVname))
        Dset = h5u.DsetCoards(OutVname, 1, [ Nz ])
        Dset.Build(OutFile, Var, Zdim)

        print("")

    InFile.close()
    OutFile.close()
