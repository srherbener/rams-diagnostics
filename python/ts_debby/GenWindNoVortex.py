#!/usr/bin/env python3
#
# Script to remove vortex from the horizontal wind fields. The method is based on that presented in
# Kurihara et al. (1993, 1995).
#   1) Separate the total wind field into a basic field and disturbance field. Use the smoothing
#      described in Kurihara et al. (1993).
#   2) Remove the vortex using a circular filter that tapers to zero at the circle edge in order
#      to prevent a discontinuity in the non-hurricane disturbance field when the vortex piece
#      is removed.
#   3) The separated field then becomes the basic (smoothed) field plus the non-hurricane
#      piece of the disturbance field.
#

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
#import scipy.interpolate as intrp
import scipy.ndimage.interpolation as intrp
import ConfigTsd as conf
import Hdf5Utils as h5u

################################
# Functions
################################

####################################################################################################
# KbrSmoothing
#
# This function will apply the smoothing described in Kurihara et al. (1993)
# applied to one field.
# 
def KbrSmoothing(Var):
    # Var is (z,y,x)

    # Use the same K coefficients as in Kurihara et al. (1993). These are in
    # terms of wavenumbers, m, which are relative to the grid spacing.
    m = np.array([ 2.0, 3.0, 4.0, 2.0, 5.0, 6.0, 7.0, 2.0, 8.0, 9.0, 2.0 ])
    K = 0.5 / ( 1.0 -  np.cos((2.0 * np.pi) / m))
    Nk = len(K)

    # First, smooth in the x direction (last dimension)
    #
    # Smoothing in KBR is on a 1 degree grid, whereas the TS Debby "lite"
    # data is on a 1/3 degree grid. Zoom down by 1/3 along the x and
    # y dimensions in order to get to approx. 1 degree grid. This way,
    # wavelengths up to 9 degrees can be filtered out.
    #
    # At the fine grid, in order to filter out 9 degree wavelengths, the
    # list in m needs to be extended to add in higher wavenumbers,
    # but doing this causes numerical issues (the K values blow up).
    #
    VarSmooth = intrp.zoom(Var,(1.0, 1.0/3.0, 1.0/3.0), order=3)
    for i in range(Nk):
       VarIm1 = VarSmooth[:,:,:-2]
       VarI   = VarSmooth[:,:,1:-1]
       VarIp1 = VarSmooth[:,:,2:]
       VarTemp = VarSmooth[:,:,1:-1] + K[i] * (VarIm1 + VarIp1 - (2.0 * VarI))
       VarSmooth[:,:,1:-1] = VarTemp

    del VarIm1
    del VarI
    del VarIp1

    # Then, smooth in the y direction (middle dimension)
    for i in range(Nk):
       VarJm1 = VarSmooth[:,:-2,:]
       VarJ   = VarSmooth[:,1:-1,:]
       VarJp1 = VarSmooth[:,2:,:]
       VarTemp = VarSmooth[:,1:-1,:] + K[i] * (VarJm1 + VarJp1 - (2.0 * VarJ))
       VarSmooth[:,1:-1,:] = VarTemp

    del VarJm1
    del VarJ
    del VarJp1

    # Return VarSmooth to the original 1/3 degree grid.
    VarTemp = np.copy(VarSmooth)
    VarSmooth = intrp.zoom(VarSmooth, (1.0,3.0,3.0), order=3)

    # Check to see if we got the original size in the resulting
    # smoothed variable. If not, then repeat the end planes
    # of the field to fill in the missing data.
    Nz, Ny, Nx = Var.shape
    Nsz, Nsy, Nsx = VarSmooth.shape
    Xdiff = Nx - Nsx
    Ydiff = Ny - Nsy

    if (Xdiff > 0):
        # VarSmooth is too small, repeat the y-z plane at the end of
        # the x-dimension.
        for i in range(Xdiff):
            # Note the use of Nz instead of Nsz. These two should match since the zoom
            # factors on the z-dimension were always set to 1.0.
            VarSmooth = np.append(VarSmooth, VarSmooth[:,:,-1].reshape((Nz,Nsy,1)), axis=2)
    elif (Xdiff < 0):
        # VarSmooth is too big, cut off the end of the x-dimension
        VarSmooth = np.copy(VarSmooth[:,:,:Xdiff])
      
    if (Ydiff > 0):
        # VarSmooth is too small, repeat the x-z plane at the end of
        # the y-dimension.
        for i in range(Ydiff):
            # Note the use of Nx instead of Nsx in the reshape() method. This is done
            # since the code above has already modified VarSmooth x-dimension to match
            # that of the input Var.
            VarSmooth = np.append(VarSmooth, VarSmooth[:,-1,:].reshape((Nz,1,Nx)), axis=1)
    elif (Ydiff < 0):
        # VarSmooth is too big, cut off the end of the x-dimension
        VarSmooth = np.copy(VarSmooth[:,:,:Ydiff])
      

    return VarSmooth


################################
# Main
################################

Tstring = conf.SetTimeString()

SimList = [
#    'TSD_SAL_DUST',
#    'TSD_SAL_NODUST',
#    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'  # vertical coords: height
Pname = '/p_coords'  # vertical coords: pressure
Tname = '/t_coords'

R0 = 5 # degrees lat/lon (~500 km)
L = R0 / 5

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Removing vortex from horizontal winds: {0:s}".format(Sim))
    print('')

    # input file and dataset names
    UinFname = "HDF5/{0:s}/HDF5/u_lite-{0:s}-AS-2006-08-20-120000-g3.h5".format(Sim, Sim)
    VinFname = "HDF5/{0:s}/HDF5/v_lite-{0:s}-AS-2006-08-20-120000-g3.h5".format(Sim, Sim)
#    UinFname = "HDF5/{0:s}/HDF5/u-{0:s}-AS-2006-08-20-120000-g3.h5".format(Sim, Sim)
#    VinFname = "HDF5/{0:s}/HDF5/v-{0:s}-AS-2006-08-20-120000-g3.h5".format(Sim, Sim)

    UinVname = '/u'
    VinVname = '/v'

    # output file and dataset names
    OutFname = "test_no_vortex_lite_{0:s}.h5".format(Sim)

    UorigVname  = '/u_orig'
    UbasicVname = '/u_basic'
    UdistVname  = '/u_dist'

    VorigVname  = '/v_orig'
    VbasicVname = '/v_basic'
    VdistVname  = '/v_dist'

    Ufile  = h5py.File(UinFname, mode='r')
    Vfile  = h5py.File(VinFname, mode='r')
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
        Ubasic = KbrSmoothing(U)
        Vbasic = KbrSmoothing(V)

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

