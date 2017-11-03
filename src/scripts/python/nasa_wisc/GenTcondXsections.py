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

# One argument which is the name of the simulation.
Sim = sys.argv[1]

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

RtFileTemplate = "SIMS/<SIM>/HDF5/rtotal-<SIM>-LS-2012-01-01-000000-g1.h5"
RtVname = "/rtotal_orig"

RvFileTemplate = "SIMS/<SIM>/HDF5/vapor-<SIM>-LS-2012-01-01-000000-g1.h5"
RvVname = "/vapor"

OutFileTemplate = "DIAGS/tcond_xsections_<SIM>.h5"
OutVnameBase = "/tcond"
LonSliceVname = "{0:s}_lon_slice".format(OutVnameBase)
LatSliceVname = "{0:s}_lat_slice".format(OutVnameBase)
LonAvgVname = "{0:s}_lon_avg".format(OutVnameBase)
LatAvgVname = "{0:s}_lat_avg".format(OutVnameBase)

# MAIN
print("************************************************************************************")
print("Generating total condensate cross sections for simulation: {0:s}".format(Sim))
print("")

# Open the input file
RtFname = RtFileTemplate.replace("<SIM>", Sim)
RtFile  = h5py.File(RtFname, mode='r')
print("  Reading: {0:s} ({1:s})".format(RtFname, RtVname))

RvFname = RvFileTemplate.replace("<SIM>", Sim)
RvFile  = h5py.File(RvFname, mode='r')
print("    Reading: {0:s} ({1:s})".format(RvFname, RvVname))

# Open the output file
OutFname = OutFileTemplate.replace("<SIM>", Sim)
OutFile = h5py.File(OutFname, mode='w')

# Read in coordinates, and transfer to the output file
print("    Reading: {0:s} ({1:s})".format(RtFname, Xname))
X = RtFile[Xname][...]
Nx = len(X)
Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
Xdim.Build(OutFile, X)

print("    Reading: {0:s} ({1:s})".format(RtFname, Yname))
Y = RtFile[Yname][...]
Ny = len(Y)
Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
Ydim.Build(OutFile, Y)

print("    Reading: {0:s} ({1:s})".format(RtFname, Zname))
Z = RtFile[Zname][...]
Nz = len(Z)
Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
Zdim.Build(OutFile, Z)

print("    Reading: {0:s} ({1:s})".format(RtFname, Tname))
T = RtFile[Tname][...]
Nt = len(T)
Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
Tdim.Build(OutFile, T)

# Create the output datasets
LonSliceDset = h5u.DsetCoards(LonSliceVname, 3, [ Nt, Nz, Ny ], chunks=(1, Nz, Ny))
LonSliceDset.Create(OutFile)

LatSliceDset = h5u.DsetCoards(LatSliceVname, 3, [ Nt, Nz, Nx ], chunks=(1, Nz, Nx))
LatSliceDset.Create(OutFile)

LonAvgDset = h5u.DsetCoards(LonAvgVname, 3, [ Nt, Nz, Ny ], chunks=(1, Nz, Ny))
LonAvgDset.Create(OutFile)

LatAvgDset = h5u.DsetCoards(LatAvgVname, 3, [ Nt, Nz, Nx ], chunks=(1, Nz, Nx))
LatAvgDset.Create(OutFile)

# Read in input variable and form the domain average
# Do one time step at a time in order to reduce memory usage
for i in range(2): # range(Nt):

    # total condensate is Rtotal - Rvapor
    # var is (z,y,x)
    TCOND = np.squeeze(RtFile[RtVname][i,...]) - np.squeeze(RvFile[RvVname][i,...])

    # Take a slice along the middle in both zonal and meridional directions
    # Ditto for averages
    MidX = int((X[0] + X[-1]) * 0.5)
    MidY = int((Y[0] + Y[-1]) * 0.5)

    TCOND_LON_SLICE = np.squeeze(TCOND[:,:,MidX])
    TCOND_LAT_SLICE = np.squeeze(TCOND[:,MidY,:])

    TCOND_LON_AVG = np.squeeze(np.mean(TCOND, axis=2))
    TCOND_LAT_AVG = np.squeeze(np.mean(TCOND, axis=1))

    # Add to output datasets
    OutFile[LonSliceVname][i,...] = TCOND_LON_SLICE
    OutFile[LatSliceVname][i,...] = TCOND_LAT_SLICE

    OutFile[LonAvgVname][i,...] = TCOND_LON_AVG
    OutFile[LatAvgVname][i,...] = TCOND_LAT_AVG

    if ((i % 10) == 0):
        print("    Processing time step: {0:d}".format(i))

print("")


RtFile.close()
RvFile.close()

# Attach dimensions to output datasets
print("  Writing {0:s} ({1:s})".format(OutFname, LonSliceVname))
print("  Writing {0:s} ({1:s})".format(OutFname, LatSliceVname))

print("  Writing {0:s} ({1:s})".format(OutFname, LonAvgVname))
print("  Writing {0:s} ({1:s})".format(OutFname, LatAvgVname))

LonSliceDset.AttachDims(OutFile, Tdim, Zdim, Ydim)
LatSliceDset.AttachDims(OutFile, Tdim, Zdim, Xdim)

LonAvgDset.AttachDims(OutFile, Tdim, Zdim, Ydim)
LatAvgDset.AttachDims(OutFile, Tdim, Zdim, Xdim)

OutFile.close()
