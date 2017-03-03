#!/usr/bin/env python3
#
# Script to remove vortex from the horizontal wind fields. This is done by smoothing the
# wind field in the region of the vortex.
#

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import scipy.ndimage.interpolation as intrp
import scipy.ndimage.filters as filt
import ConfigTsd as conf
import Hdf5Utils as h5u

################################
# Functions
################################

####################################################################################################
# KbrSmoothing
#
# This function applies the vortex smoothing on the horizontal wind field using the method
# outlined in Kurihara et al. (1993).
# 
def KbrSmoothing(Var, VarSelect):
    # Use the same K coefficients as in Kurihara et al. (1993). These are in
    # terms of wavenumbers, m, which are relative to the grid spacing.
    m = np.array([ 2.0, 3.0, 4.0, 2.0, 5.0, 6.0, 7.0, 2.0, 8.0, 9.0, 2.0 ])
    K = 0.5 / ( 1.0 -  np.cos((2.0 * np.pi) / m))
    Nk = len(K)

    # First, smooth in the x direction (last dimension)
    VarSmooth = np.copy(Var)
    for i in range(Nk):
        VarTemp = np.copy(VarSmooth)

        VarIm1 = VarTemp[:,:,:-2]
        VarI   = VarTemp[:,:,1:-1]
        VarIp1 = VarTemp[:,:,2:]
        VarSmooth[:,:,1:-1] = VarI + K[i] * (VarIm1 + VarIp1 - (2.0 * VarI))
        # Take difference and only apply difference where VarSelect has 1.0 values
        # (VarSelect has either 0.0 or 1.0 values.)
        VarDiff = (VarSmooth - VarTemp) * VarSelect
        VarSmooth = VarTemp + VarDiff
        

    del VarIm1
    del VarI
    del VarIp1
    
    # Then, smooth in the y direction (middle dimension)
    for i in range(Nk):
        VarTemp = np.copy(VarSmooth)

        VarJm1 = VarTemp[:,:-2,:]
        VarJ   = VarTemp[:,1:-1,:]
        VarJp1 = VarTemp[:,2:,:]
        VarSmooth[:,1:-1,:] = VarJ + K[i] * (VarJm1 + VarJp1 - (2.0 * VarJ))
        # Take difference and only apply difference where VarSelect has 1.0 values
        # (VarSelect has either 0.0 or 1.0 values.)
        VarDiff = (VarSmooth - VarTemp) * VarSelect
        VarSmooth = VarTemp + VarDiff

    del VarJm1
    del VarJ
    del VarJp1
    
    return VarSmooth

####################################################################################################
# UniformSmoothing
#
def UniformSmoothing(Var, VarSelect):

    Fsize = 5

    VarSmooth = np.copy(Var)
    for i in range(10):
        VarTemp = np.copy(VarSmooth)
        VarSmooth = filt.uniform_filter(VarTemp, size=Fsize, mode='constant')
        VarDiff = (VarSmooth - VarTemp) * VarSelect
        VarSmooth = VarTemp + VarDiff

    return VarSmooth


####################################################################################################
# VortexSmoothing
#
# This function applies the vortex smoothing on the horizontal wind field.
#
# The 'method' argument selects the smoothing algortime.
#   'kbr': smoothing as described in Kurihara et al. (1993)
#   'uni': uniform 2D filter from scipy
# 
def VortexSmoothing(X, Y, Var, Method, MaxRadius, StormLocX, StormLocY):
    # Var is organized as (z,y,x), and we are smoothing regions in the (x,y) planes
    # Construct an array with same shape as Var that contains the radial lengths
    # from the storm center. Use this to limit the smoothing to just the storm
    # region.
    Nz, Ny, Nx = Var.shape
    StormX = X - StormLocX
    StormY = Y - StormLocY
    [ Xgrid, Ygrid ] = np.meshgrid(StormX, StormY)
    Radius = np.sqrt(np.square(Xgrid) + np.square(Ygrid))
    Radius = np.tile(Radius, (Nz,1,1))

    # Create an array with zeros outside the storm, and ones inside the storm
    VarSelect = np.zeros([ Nz, Ny, Nx])
    VarSelect[Radius <= MaxRadius] = 1.0

    if (Method == 'kbr'):
        VarSmooth = KbrSmoothing(Var, VarSelect)
    elif (Method == 'uni'):
        VarSmooth = UniformSmoothing(Var, VarSelect)

    return VarSmooth


################################
# Main
################################

Tstring = conf.SetTimeString()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'  # vertical coords: height
Pname = '/p_coords'  # vertical coords: pressure
Tname = '/t_coords'

ScriptName = sys.argv[0]
Usage = "USAGE: {0:s} <smooth_method> <vert_coord_type>".format(ScriptName)
if (len(sys.argv) != 3):
    print("ERROR: {0:s}: Must supply two arguments".format(ScriptName))
    print(Usage)
    sys.exit(2)

SmoothMethod = sys.argv[1]
VcoordType   = sys.argv[2]

if ((SmoothMethod != 'uni') and (SmoothMethod != 'kbr')):
    print("ERROR: {0:s}: <smooth_method>, ({1:s}), must be one of 'uni', 'kbr'".format(ScriptName, SmoothMethod))
    print(Usage)
    sys.exit(3)

if ((VcoordType != 'z') and (VcoordType != 'p')):
    print("ERROR: {0:s}: <vert_coord_type>, ({1:s}), must be one of 'z', 'p'".format(ScriptName, VcoordType))
    print(Usage)
    sys.exit(4)

if (VcoordType == 'z'):
    FileSuffix = "-AS-2006-08-20-120000-g3.h5"
elif (VcoordType == 'p'):
    FileSuffix = "-AP-2006-08-20-120000-g3.h5"

ScFileSuffix = "-AS-2006-08-20-120000-g3.h5" # this only depends on horizontal dimensions

MaxRadius = 3 # degrees (roughly 300 km)

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Removing vortex from horizontal winds: {0:s}".format(Sim))
    print("  Vortex smoothing method: {0:s}".format(SmoothMethod))
    print("  Vertical coordinate type: {0:s}".format(VcoordType))
    print("  Maximum radius: {0:.1f}".format(MaxRadius))
    print('')

    # input file and dataset names
    UinFname = "HDF5/{0:s}/HDF5/u_lite-{0:s}{1:s}".format(Sim, FileSuffix)
    VinFname = "HDF5/{0:s}/HDF5/v_lite-{0:s}{1:s}".format(Sim, FileSuffix)
    ScFname  = "HDF5/{0:s}/HDF5/storm_center-{0:s}{1:s}".format(Sim, ScFileSuffix)

    UinVname = '/u'
    VinVname = '/v'
    SxlocVname = '/press_cent_xloc'
    SylocVname = '/press_cent_yloc'

    # output file and dataset names
    OutFname = "HDF5/{0:s}/HDF5/hwinds_no_vortex_lite-{0:s}{1:s}".format(Sim, FileSuffix)

    UorigVname  = '/u_orig'
    UbasicVname = '/u_basic'
    UdistVname  = '/u_dist'

    VorigVname  = '/v_orig'
    VbasicVname = '/v_basic'
    VdistVname  = '/v_dist'

    Ufile  = h5py.File(UinFname, mode='r')
    Vfile  = h5py.File(VinFname, mode='r')
    ScFile = h5py.File(ScFname, mode='r')
    OutFile = h5py.File(OutFname, mode='w')

    # Build the coordinates
    print("  Reading: {0:s} ({1:s})".format(UinFname, Xname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, Xname))
    X = Ufile[Xname][...]
    Nx = len(X)
    Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
    Xdim.Build(OutFile, X)

    print("  Reading: {0:s} ({1:s})".format(UinFname, Yname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, Yname))
    Y = Ufile[Yname][...]
    Ny = len(Y)
    Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
    Ydim.Build(OutFile, Y)

    print("  Reading: {0:s} ({1:s})".format(UinFname, Zname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, Zname))
    Z = Ufile[Zname][...]
    Nz = len(Z)
    Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
    Zdim.Build(OutFile, Z)

    print("  Reading: {0:s} ({1:s})".format(UinFname, Tname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, Tname))
    T = Ufile[Tname][...]
    Nt = len(T)
    Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
    Tdim.Build(OutFile, T)

    # Read the storm center values
    print("  Reading: {0:s} ({1:s})".format(ScFname, SxlocVname))
    print("  Reading: {0:s} ({1:s})".format(ScFname, SylocVname))
    StormLocX = ScFile[SxlocVname][...]
    StormLocY = ScFile[SylocVname][...]

    print("")

    # Create the output datasets
    UorigDset = h5u.DsetCoards(UorigVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
    VorigDset = h5u.DsetCoards(VorigVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
    UbasicDset = h5u.DsetCoards(UbasicVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
    VbasicDset = h5u.DsetCoards(VbasicVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
    UdistDset = h5u.DsetCoards(UdistVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
    VdistDset = h5u.DsetCoards(VdistVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))

    UorigDset.Create(OutFile)
    VorigDset.Create(OutFile)
    UbasicDset.Create(OutFile)
    VbasicDset.Create(OutFile)
    UdistDset.Create(OutFile)
    VdistDset.Create(OutFile)

    print("  Reading: {0:s} ({1:s})".format(UinFname, UinVname))
    print("  Reading: {0:s} ({1:s})".format(VinFname, VinVname))
    print("")

    # In order to keep memory requirements down, process one time step at a time
    for it in range(Nt):

        # U and V are dimensioned as (t,z,y,x). So slice on the first dimension.
        U = Ufile[UinVname][it,:,:,:]
        V = Vfile[VinVname][it,:,:,:]
        
        # Basic fields (smoothed)
        Ubasic = VortexSmoothing(X, Y, U, SmoothMethod, MaxRadius, StormLocX[it], StormLocY[it])
        Vbasic = VortexSmoothing(X, Y, V, SmoothMethod, MaxRadius, StormLocX[it], StormLocY[it])

        # Disturbance field is total - basic fields
        Udist = U - Ubasic
        Vdist = V - Vbasic

        # Write fields into output file
        OutFile[UorigVname][it,:,:,:] = U
        OutFile[VorigVname][it,:,:,:] = V

        OutFile[UbasicVname][it,:,:,:] = Ubasic
        OutFile[VbasicVname][it,:,:,:] = Vbasic

        OutFile[UdistVname][it,:,:,:] = Udist
        OutFile[VdistVname][it,:,:,:] = Vdist
        
        
        if ((it % 10) == 0):
           print("  Processing time step: {0:d}".format(it))

    print("")

    # Attach the dimensions to U and V.
    print("  Writing: {0:s} ({1:s})".format(OutFname, UorigVname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, VorigVname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, UbasicVname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, VbasicVname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, UdistVname))
    print("  Writing: {0:s} ({1:s})".format(OutFname, VdistVname))
    UorigDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)
    VorigDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)
    UbasicDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)
    VbasicDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)
    UdistDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)
    VdistDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)

    Ufile.close()
    Vfile.close()
    OutFile.close()

    print("")

