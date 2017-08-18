#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))

import h5py
import glob
import numpy as np
import ConfigRce as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()

Lv = 2.5e6
Cp = 1004

# RAMS files are lite version that do not contain the coordinates, so use the functions
# from the config utilities.
X = conf.SetHcoords(-2.15497, 2.15497, 480)
Nx = len(X)

Y = conf.SetHcoords(-2.15497, 2.15497, 480)
Ny = len(Y)

Zall = conf.SetZcoords()
NzAll = len(Zall)
DeltaZ = []
for iz in range(NzAll-1):
    DeltaZ.append(Zall[iz+1] - Zall[iz])

Xname = "/x_coords"
Yname = "/y_coords"
Zname = "/z_coords"
Tname = "/t_coords"


SimList = [
    'RCE_1km',
    'RCE_1km_SM',
    'RCE_1km_DM',
    'RCE_1km_DP',
    ]
Nsims = len(SimList)

VarList = [
    'pcp_water',
    'olr',
    'pcp_rate',
    'sfc_lhf',
    'sfc_shf',
    ]
Nvars = len(VarList)

### TEST #### InFilePatternTemplate = "SIMDATA/<SIM>/RAMS/<SIM>-L-2012-01-01-0[0-2]*.h5"
InFilePatternTemplate = "SIMDATA/<SIM>/RAMS/<SIM>-L-2012-*.h5"
OutFileTemplate = "SIMDATA/<SIM>/HDF5/<VNAME>-<SIM>-AC-2012-01-01-000000-g1.h5"

for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Extracting RCE variables for simulation: {0:s}".format(Sim))
    print("")

    # Input file
    InFilePattern = InFilePatternTemplate.replace("<SIM>", Sim)
    FileList = glob.glob(InFilePattern)
    Nt = len(FileList)
    T = np.arange(Nt) * 3600.0 # sim time in seconds, start at zero

    # Open all of the output files
    OutFnames = []
    OutFiles = []
    for ivar in range(Nvars):
        Vname = VarList[ivar]

        OutFnames.append(OutFileTemplate.replace("<SIM>", Sim).replace("<VNAME>", Vname))
        OutFiles.append(h5py.File(OutFnames[ivar], mode='w'))

    # Process one input file at a time, read in all vars and write into
    # output files.
    it = 0
    OutDsets = []
    for InFname in sorted(FileList):
        print("  Reading: {0:s}, Time step: {1:d}".format(InFname, it+1))
        InFile = h5py.File(InFname, mode='r')

        for ivar in range(Nvars):
            Vname = VarList[ivar]
            print("    Variable: {0:s}".format(Vname))

            if (Vname == 'pcp_water'):
                DN0 = np.squeeze(InFile['DN0'][...])  # kg/m^3
                RV = np.squeeze(InFile['RV'][...])    # kg/kg
                Q = RV * DN0
                VAR = np.zeros((Ny,Nx))   # mm
                for iz in range(NzAll-1):
                    VAR = VAR + (DeltaZ[iz] * Q[iz+1,...])

            elif (Vname == 'olr'):
                VAR = np.squeeze(InFile['LWUP'][-1,...])  # W/m^2

            elif (Vname == 'pcp_rate'):
                VAR = np.squeeze(InFile['PCPRR'][...]) * 86400  # mm/day (86400 seconds in a day)

            elif (Vname == 'sfc_lhf'):
                VAR = np.squeeze(InFile['SFLUX_R'][...]) * Lv  # W/m^2

            elif (Vname == 'sfc_shf'):
                VAR = np.squeeze(InFile['SFLUX_T'][...]) * Cp  # W/m^2

            else:
                print("      Warning: undefined variable ({0:s}), skipping this variable".format(Vname))
                continue

            # Cast a 2D variable into a 3D variable with a z-dimension length of 1
            if (len(VAR.shape) == 2):
                Z = [ 0 ]
                Nz = 1
                VAR = VAR.reshape((1, VAR.shape[0], VAR.shape[1]))
            else:
                Z = Zall
                Nz = NzAll

            # If this is the first time step, build the coordinates and dataset for the variable
            if (it == 0):
                # Build the coordinates in the output file
                Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
                Xdim.Build(OutFiles[ivar], X)

                Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
                Ydim.Build(OutFiles[ivar], Y)

                Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
                Zdim.Build(OutFiles[ivar], Z)

                Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
                Tdim.Build(OutFiles[ivar], T)

                OutDsets.append(h5u.DsetCoards(Vname, 4, (Nt, Nz, Ny, Nx), chunks=(1, Nz, Ny, Nx)))
                OutDsets[ivar].Create(OutFiles[ivar])

            # Write the var into the output file
            OutFiles[ivar][Vname][it,...] = VAR

        InFile.close()
        it = it + 1

    print("")

    # Attach dims and close the output files
    for ivar in range(Nvars):
        print("  Writing: {0:s} ({1:s})".format(OutFnames[ivar], VarList[ivar]))
        OutDsets[ivar].AttachDims(OutFiles[ivar], Tdim, Zdim, Ydim, Xdim)
        OutFiles[ivar].close()

    print("")

