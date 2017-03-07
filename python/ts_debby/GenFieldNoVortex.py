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
import scipy.signal as sig
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
# ShumanSmoothing
#
# This function will do 2D smoothing on horizontal layers according to the technique described
# in Shuman (1957).
#
# Smoothing is done in 3 stages using:
#     Nu = 0.49965
#     Nu = -0.22227 + 0.64240j
#     Nu = -0.22227 - 0.64240j
#
# Formula for each stage is:
#    Zbar = Z0 + 1/2*Nu*(1-Nu)*(Z2+Z4+Z6+Z8 - 4*Z0) + 1/4*Nu**2*(Z1+Z3+Z5+Z7 - 4*Z0)
#
#    Where the locations of the Z's are (on the horizontal plane):
#
#       Z1  Z8  Z7
#       Z2  Z0  Z6
#       Z3  Z4  Z5
#
def NuConvolve(Var, Nu):
    Nz, Ny, Nx = Var.shape

    NuCoeffEven = 0.5 * Nu * (1.0-Nu)
    NuCoeffOdd  = 0.25 * Nu * Nu
    NuCoeffCntr = 1 - (4.0 * (NuCoeffEven + NuCoeffOdd))

    Coeffs = np.array([
      [ NuCoeffOdd,  NuCoeffEven, NuCoeffOdd  ],
      [ NuCoeffEven, NuCoeffCntr, NuCoeffEven ],
      [ NuCoeffOdd,  NuCoeffEven, NuCoeffOdd  ] ])

    VarSmooth = np.zeros([Nz, Ny, Nx], dtype=np.complex_)
    for k in range(Nz):
        VarSmooth[k,:,:] = sig.convolve2d(np.squeeze(Var[k,:,:]), Coeffs, mode='same', boundary='symm')

    return VarSmooth


def ShumanSmoothing(Var, VarSelect):
    # Utilize a 9-point mesh (3x3 array with current point in middle) for the smoothing. 

    Nz, Ny, Nx = Var.shape

    a = 0.49965
    b = -0.22227+0.64240j
    c = np.conj(b)
    Nu = np.array([ a, b, c ])

    VarSmooth = np.copy(Var)
    for n in range(10):
        VarTemp = np.copy(VarSmooth)
        # do smoothing across the horizontal
        for k in range(Nz):
            VarSmooth = NuConvolve(VarTemp, Nu[0])
            VarSmooth = NuConvolve(VarSmooth, Nu[1])
            VarSmooth = NuConvolve(VarSmooth, Nu[2])
        # VarSmooth is type complex now. The set of Nu's are required to contain complex
        # conjugates which means that VarSmooth is supposed to be real (imaginary part
        # equal to zero) at the end of the smoothing. However, discrete numerics will
        # likely leave tiny non-zero imaginary parts. Take the real part of VarSmooth
        # in the statement below to throw away the near zero imaginary parts.
        VarDiff = (VarSmooth.real - VarTemp) * VarSelect
        VarSmooth = VarTemp + VarDiff

    return VarSmooth


####################################################################################################
# UniformSmoothing
#
def UniformSmoothing(Var, VarSelect):

    Fsize = 5

    Nz, Ny, Nx = Var.shape

    VarSmooth = np.copy(Var)
    for i in range(10):
        VarTemp = np.copy(VarSmooth)
        # do the smoothing only across the horizontal
        for k in range(Nz):
            VarSmooth[k,:,:] = filt.uniform_filter(np.squeeze(VarTemp[k,:,:]), size=Fsize, mode='reflect')
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
    elif (Method == 'shu'):
        VarSmooth = ShumanSmoothing(Var, VarSelect)
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

VarList = [
    [ "HDF5/<SIM>/HDF5/u_lite-<SIM>-AS-2006-08-20-120000-g3.h5", '/u', "HDF5/<SIM>/HDF5/u_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "/press_cent_xloc", "/press_cent_yloc" ],
    [ "HDF5/<SIM>/HDF5/v_lite-<SIM>-AS-2006-08-20-120000-g3.h5", '/v', "HDF5/<SIM>/HDF5/v_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "/press_cent_xloc", "/press_cent_yloc" ],

    [ "HDF5/<SIM>/HDF5/tempc_lite-<SIM>-AS-2006-08-20-120000-g3.h5", '/tempc', "HDF5/<SIM>/HDF5/tempc_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "/press_cent_xloc", "/press_cent_yloc" ],

    [ "HDF5/<SIM>/HDF5/vapor_lite-<SIM>-AS-2006-08-20-120000-g3.h5", '/vapor', "HDF5/<SIM>/HDF5/vapor_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "/press_cent_xloc", "/press_cent_yloc" ],

    [ "HDF5/<SIM>/HDF5/u_lite-<SIM>-AP-2006-08-20-120000-g3.h5", '/u', "HDF5/<SIM>/HDF5/u_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5", "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "/press_cent_xloc", "/press_cent_yloc" ],
    [ "HDF5/<SIM>/HDF5/v_lite-<SIM>-AP-2006-08-20-120000-g3.h5", '/v', "HDF5/<SIM>/HDF5/v_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5", "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "/press_cent_xloc", "/press_cent_yloc" ],

    [ "HDF5/<SIM>/HDF5/tempc_lite-<SIM>-AP-2006-08-20-120000-g3.h5", '/tempc', "HDF5/<SIM>/HDF5/tempc_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5", "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "/press_cent_xloc", "/press_cent_yloc" ],

    [ "HDF5/<SIM>/HDF5/vapor_lite-<SIM>-AP-2006-08-20-120000-g3.h5", '/vapor', "HDF5/<SIM>/HDF5/vapor_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5", "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5", "/press_cent_xloc", "/press_cent_yloc" ],
    ]
Nvars = len(VarList)

Xname = "/x_coords"
Yname = "/y_coords"
Zname = "/z_coords"
Tname = "/t_coords"

MaxRadius = 3 # degrees (roughly 300 km)
# Shuman method doesn't seem to filter enough
SmoothMethod = "uni"

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Removing vortex from field: {0:s}".format(Sim))
    print("  Vortex smoothing method: {0:s}".format(SmoothMethod))
    print("  Maximum radius: {0:.1f}".format(MaxRadius))
    print("")

    for ivar in range(Nvars):
        InFname      = VarList[ivar][0].replace("<SIM>", Sim)
        InVname      = VarList[ivar][1]
        OutFname     = VarList[ivar][2].replace("<SIM>", Sim)
        ScFname      = VarList[ivar][3].replace("<SIM>", Sim)
        ScXlocVname  = VarList[ivar][4]
        ScYlocVname  = VarList[ivar][5]
    
        # output dataset names
        OrigVname  = "{0:s}_orig".format(InVname)
        BasicVname = "{0:s}_basic".format(InVname)
        DistVname  = "{0:s}_dist".format(InVname)
    
        InFile  = h5py.File(InFname, mode='r')
        ScFile = h5py.File(ScFname, mode='r')
        OutFile = h5py.File(OutFname, mode='w')
    
        # Build the coordinates
        print("  Reading: {0:s} ({1:s})".format(InFname, Xname))
        print("  Writing: {0:s} ({1:s})".format(OutFname, Xname))
        X = InFile[Xname][...]
        Nx = len(X)
        Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
        Xdim.Build(OutFile, X)
    
        print("  Reading: {0:s} ({1:s})".format(InFname, Yname))
        print("  Writing: {0:s} ({1:s})".format(OutFname, Yname))
        Y = InFile[Yname][...]
        Ny = len(Y)
        Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
        Ydim.Build(OutFile, Y)
    
        print("  Reading: {0:s} ({1:s})".format(InFname, Zname))
        print("  Writing: {0:s} ({1:s})".format(OutFname, Zname))
        Z = InFile[Zname][...]
        Nz = len(Z)
        Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
        Zdim.Build(OutFile, Z)
    
        print("  Reading: {0:s} ({1:s})".format(InFname, Tname))
        print("  Writing: {0:s} ({1:s})".format(OutFname, Tname))
        T = InFile[Tname][...]
        Nt = len(T)
        Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
        Tdim.Build(OutFile, T)
    
        # Read the storm center values
        print("  Reading: {0:s} ({1:s})".format(ScFname, ScXlocVname))
        print("  Reading: {0:s} ({1:s})".format(ScFname, ScYlocVname))
        StormLocX = ScFile[ScXlocVname][...]
        StormLocY = ScFile[ScYlocVname][...]
    
        print("")
    
        # Create the output datasets
        OrigDset = h5u.DsetCoards(OrigVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
        BasicDset = h5u.DsetCoards(BasicVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
        DistDset = h5u.DsetCoards(DistVname, 4, [ Nt, Nz, Ny, Nx ], chunks=( 1, Nz, Ny, Nx ))
    
        OrigDset.Create(OutFile)
        BasicDset.Create(OutFile)
        DistDset.Create(OutFile)
    
        print("")
    
        # In order to keep memory requirements down, process one time step at a time
        print("  Reading: {0:s} ({1:s})".format(InFname, InVname))
        for it in range(Nt):
    
            # Input field is dimensioned as (t,z,y,x). So slice on the first dimension.
            Var = InFile[InVname][it,:,:,:]
            
            # Basic field (smoothed)
            VarBasic = VortexSmoothing(X, Y, Var, SmoothMethod, MaxRadius, StormLocX[it], StormLocY[it])
    
            # Disturbance field is total - basic fields
            VarDist = Var - VarBasic
    
            # Write fields into output file
            OutFile[OrigVname][it,:,:,:] = Var
            OutFile[BasicVname][it,:,:,:] = VarBasic
            OutFile[DistVname][it,:,:,:] = VarDist
            
            
            if ((it % 10) == 0):
               print("  Processing time step: {0:d}".format(it))
    
        print("")
    
        # Attach the dimensions to output fields
        print("  Writing: {0:s} ({1:s})".format(OutFname, OrigVname))
        print("  Writing: {0:s} ({1:s})".format(OutFname, BasicVname))
        print("  Writing: {0:s} ({1:s})".format(OutFname, DistVname))
        OrigDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)
        BasicDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)
        DistDset.AttachDims(OutFile, Tdim, Zdim, Ydim, Xdim)
    
        InFile.close()
        ScFile.close()
        OutFile.close()
    
        print("")

