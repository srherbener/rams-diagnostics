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
    'RCE_3km_2mom',
    'RCE_3km_2mom_db',
    'RCE_3km_2mom_dm',
    ]
Nsims = len(SimList)

# Start and end times are time step indexes. Each time step is on hour and
# we want the final 10 days.
#   T1 = -240   240 timesteps from the end (10 X 24 hours)
#   T2 = -1     end
VarList = [
    [ "/avg_theta",      -240, -1,   1.0, 0.0 ],
    [ "/avg_total_cond", -240, -1, 1.0e3, 0.0 ],
    [ "/avg_vapor",      -240, -1, 1.0e3, 0.0 ],
    ] 
Nvars = len(VarList)

Z = conf.SetZcoords()
Nz = len(Z)
Zname = '/z_coords'

InFileTemplate = "DIAGS/dom_averages_<SIM>.h5"

# MAIN
for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating profiles for simulation: {0:s}".format(Sim))
    print("")

    # Open the output file
    OutFname = "DIAGS/profiles_{0:s}.h5".format(Sim)
    OutFile = h5py.File(OutFname, mode='w')

    # Open the input file
    InFname = InFileTemplate.replace("<SIM>", Sim)
    InFile  = h5py.File(InFname, mode='r')

    # Read in Z coordinates and build the Z dimension in the output file
    Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
    Zdim.Build(OutFile, Z)

    # Read input variables, calculate horizontal averages and write out results
    for ivar in range(Nvars):
        Vname  = VarList[ivar][0]
        T1     = VarList[ivar][1]
        T2     = VarList[ivar][2]
        Scale  = VarList[ivar][3]
        Offset = VarList[ivar][4]

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
