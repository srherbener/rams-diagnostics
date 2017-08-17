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

# RAMS files are lite version that do not contain the coordinates, so use the functions
# from the config utilities.
X = conf.SetHcoords(-2.15497, 2.15497, 480)
Nx = len(X)

Y = conf.SetHcoords(-2.15497, 2.15497, 480)
Ny = len(Y)

Z = conf.SetZcoords()
Nz = len(Z)

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

#InFilePatternTemplate = "RAMS/<SIM>/RAMS/<SIM>-L-2012-*.h5"
InFilePatternTemplate = "SIMDATA/<SIM>/RAMS/<SIM>-L-2012-01-01-0*.h5"
OutFileTemplate = "SIMDATA/<SIM>/HDF5/rce_vars-<SIM>-AC-2012-01-01-000000-g1.h5"

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

    # Open the output file
    OutFname = OutFileTemplate.replace("<SIM>", Sim)
    print("OutFname", OutFname)
    OutFile = h5py.File(OutFname, mode='w')

    # Build the coordinates in the output file
    Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
    Xdim.Build(OutFile, X)

    Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
    Ydim.Build(OutFile, Y)

    Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
    Zdim.Build(OutFile, Z)

    Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
    Tdim.Build(OutFile, T)

    


    print("")
    OutFile.close()



### 
###     # Read in input variable and form the domain average
###     # Do one time step at a time in order to reduce memory usage
###     for i in range(Nt):
###         # Read in the var and determine if 2D (t,y,x) or 3D (t,z,y,x) field.
###         Var = np.squeeze(InFile[InVname][i,...])
###         Ndims = Var.ndim
### 
###         if (Ndims == 2):
###             # Var is (y,x)
###             if (i == 0):
###                 VarAvg = np.zeros(Nt)
###                 VarAvg[i] = CalcHorizAvg(Var, Select)
###         elif (Ndims == 3):
###             # Var is (z,y,x)
###             if (i == 0):
###                 VarAvg = np.zeros((Nt,Nz))
###             for k in range(Nz):
###                 VarAvg[i,k] = CalcHorizAvg(np.squeeze(Var[i,:,:]), Select)
### 
###         if ((i % 10) == 0):
###             print("    Processing time step: {0:d}".format(i))
### 
###     print("")
### 
###     InFile.close()
### 
###     # Write out average
###     print("  Writing: {0:s} ({1:s})".format(OutFname, OutVname))
###     Dset = h5u.DsetCoards(OutVname, Ndims, VarAvg.shape)
###     if (Ndims == 2):
###         Dset.Build(OutFile, VarAvg, Tdim)
###     elif (Ndims == 3):
###         Dset.Build(OutFile, VarAvg, Tdim, Zdim)
### 
###     print("")
