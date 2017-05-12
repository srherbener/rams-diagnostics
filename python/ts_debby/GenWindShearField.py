#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

# For now, need to run either the set with full horizontal data, or the set with reduced
# horizontal data. Can't run both at the same time.
FilterFnameTemplate = "FILTERS/all_500_lite_<SIM>.h5"
FilterVname = "/filter"

UfnameTemplate = "HDF5/<SIM>/HDF5/u_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
VfnameTemplate = "HDF5/<SIM>/HDF5/v_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
Uvname = "/u_basic"
Vvname = "/v_basic"

OfnameTemplate = "DIAGS/wind_shear_nv_lite_<SIM>.h5"
UshearVname = "/u_shear"
VshearVname = "/v_shear"
MagShearVname = "/mag_shear"
AngShearVname = "/angle_shear"
MagShearStormVname = "/mag_shear_storm"
AvgMagShearStormVname = "/avg_mag_shear_storm"

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'  # vertical coords: height
Pname = '/p_coords'  # vertical coords: pressure
Tname = '/t_coords'

# For selection of levels, SFC is layer near surface
#              Press coords    Height coords
#  SAL bot        700 mb         3500 m
#  SAL top        500 mb         5500 m
#
#  SFC bot       1000 mb          110 m
#  SFC top        925 mb          761 m

ZsalBot = 700 # mb
ZsalTop = 700 # mb

ZsfcBot = 925 # mb
ZsfcTop = 925 # mb

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Creating wind shear field for simulation: {0:s}".format(Sim))
    print("")

    FilterFname = FilterFnameTemplate.replace("<SIM>", Sim)
    Ufname      = UfnameTemplate.replace("<SIM>", Sim)
    Vfname      = VfnameTemplate.replace("<SIM>", Sim)
    Ofname      = OfnameTemplate.replace("<SIM>", Sim)

    print("  Reading {0:s} ({1:s})".format(FilterFname, FilterVname))
    print("  Reading {0:s} ({1:s})".format(Vfname, Vvname))
    print("  Reading {0:s} ({1:s})".format(Vfname, Vvname))
    print("")

    # Read input data. If this is the first set, read in the coordinates
    # and build the dimensions.
    Ffile = h5py.File(FilterFname, mode='r')
    Ufile = h5py.File(Ufname, mode='r')
    Vfile = h5py.File(Vfname, mode='r')

    Ofile = h5py.File(Ofname, mode='w')

    print("    Reading {0:s} ({1:s})".format(Ufname, Xname))
    X = Ufile[Xname][...]
    Nx = len(X)
    Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
    Xdim.Build(Ofile, X)

    print("    Reading {0:s} ({1:s})".format(Ufname, Yname))
    Y = Ufile[Yname][...]
    Ny = len(Y)
    Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
    Ydim.Build(Ofile, Y)

    print("    Reading {0:s} ({1:s})".format(Ufname, Zname))
    Z = Ufile[Zname][...]
    Nz = len(Z)
    Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
    Zdim.Build(Ofile, Z)

    print("    Reading {0:s} ({1:s})".format(Ufname, Tname))
    T = Ufile[Tname][...]
    Nt = len(T)
    Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
    Tdim.Build(Ofile, T)

    print("")

    # Take the layer averages of the SAL level winds, and near surface winds,
    # and use these averages to caclulate the shear (vector difference).

    # Find the indices for the top and bottom of the layers.
    # pressure coords 
    Select = np.where((Z <= ZsfcBot) * (Z >= ZsfcTop))
    SFC_Z1 = Select[0][0]
    SFC_Z2 = Select[0][-1]

    Select = np.where((Z <= ZsalBot) * (Z >= ZsalTop))
    SAL_Z1 = Select[0][0]
    SAL_Z2 = Select[0][-1]

    # Create the output datasets using the COARDS convention so that grads
    # can read these datasets directly (sdfopen).
    UshearDset = h5u.DsetCoards(UshearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    VshearDset = h5u.DsetCoards(VshearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    MagShearDset = h5u.DsetCoards(MagShearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    AngShearDset = h5u.DsetCoards(AngShearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    MagShearStormDset = h5u.DsetCoards(MagShearStormVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    AvgMagShearStormDset = h5u.DsetCoards(AvgMagShearStormVname, 1, [ Nt ])

    UshearDset.Create(Ofile)
    VshearDset.Create(Ofile)
    MagShearDset.Create(Ofile)
    AngShearDset.Create(Ofile)
    MagShearStormDset.Create(Ofile)
    AvgMagShearStormDset.Create(Ofile)

    # Process one time step at a time in order to mitigate large memory allocation
    for it in range(Nt):
        FILTER = np.squeeze(Ffile[FilterVname][it,...])
        U = np.squeeze(Ufile[Uvname][it,...])
        V = np.squeeze(Vfile[Vvname][it,...])

        # Python indexing Z1:Z2 stops one before Z2. We want Z2 included so  use
        # indexing Z1:Z2+1.
        #
        # Take mean on z-dimension (1st dimension) which will yield a layer average.
        #
        # Then take the difference between layers to get the u,v components of the shear vectors.
        U_SFC = np.squeeze(np.mean(U[SFC_Z1:SFC_Z2+1,...],0))
        U_SAL = np.squeeze(np.mean(U[SAL_Z1:SAL_Z2+1,...],0))

        V_SFC = np.squeeze(np.mean(V[SFC_Z1:SFC_Z2+1,...],0))
        V_SAL = np.squeeze(np.mean(V[SAL_Z1:SAL_Z2+1,...],0))

        # Cartesian
        U_SHEAR = U_SAL - U_SFC
        V_SHEAR = V_SAL - V_SFC

        # Magnitude, angle
        # Create a storm relative magnitude shear field as well - this may be
        # helpful for doing comparisons. Filter has 1's in the storm region and
        # 0's otherwise.
        MAG_SHEAR = np.sqrt(np.square(U_SHEAR) + np.square(V_SHEAR))
        ANG_SHEAR = np.arctan2(V_SHEAR, U_SHEAR)
        MAG_SHEAR_STORM = np.multiply(MAG_SHEAR, FILTER)

        # Get an average shear magnitude over the storm reagion. 
        AVG_MAG_SHEAR_STORM = np.mean(MAG_SHEAR_STORM[MAG_SHEAR_STORM > 0.0])

        # Write fields into output file
        Ofile[UshearVname][it,:,:] = U_SHEAR
        Ofile[VshearVname][it,:,:] = V_SHEAR
        Ofile[MagShearVname][it,:,:] = MAG_SHEAR
        Ofile[AngShearVname][it,:,:] = ANG_SHEAR
        Ofile[MagShearStormVname][it,:,:] = MAG_SHEAR_STORM
        Ofile[AvgMagShearStormVname][it] = AVG_MAG_SHEAR_STORM

    Ufile.close()
    Vfile.close()
   
    # Attach dimensions to output fields
    print("  Writing {0:s} ({1:s})".format(Ofname, UshearVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, VshearVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShearVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, AngShearVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShearStormVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, AvgMagShearStormVname))

    UshearDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    VshearDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    MagShearDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    AngShearDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    MagShearStormDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    AvgMagShearStormDset.AttachDims(Ofile, Tdim)

    Ofile.close()
    print("")

