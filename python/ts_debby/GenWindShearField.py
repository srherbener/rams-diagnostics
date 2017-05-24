#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u

Vcoord = 'p' # set to 'z' (height) or 'p' (pressure)

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

if (Vcoord == 'p'):
  UfnameTemplate = "HDF5/<SIM>/HDF5/u_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
  VfnameTemplate = "HDF5/<SIM>/HDF5/v_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
  #UfnameTemplate = "HDF5/<SIM>/HDF5/u_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
  #VfnameTemplate = "HDF5/<SIM>/HDF5/v_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
elif (Vcoord == 'z'):
  UfnameTemplate = "HDF5/<SIM>/HDF5/u_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
  VfnameTemplate = "HDF5/<SIM>/HDF5/v_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
  #UfnameTemplate = "HDF5/<SIM>/HDF5/u_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
  #VfnameTemplate = "HDF5/<SIM>/HDF5/v_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
else:
  print("ERROR: must set Vcoord to either 'z' or 'p'")
  exit(1)

Uvname = "/u_basic"
Vvname = "/v_basic"
#Uvname = "/u"
#Vvname = "/v"

OfnameTemplate = "DIAGS/wind_shear_nv_lite_<SIM>.h5"
UshearVname = "/u_shear"
VshearVname = "/v_shear"
MagShearVname = "/mag_shear"
AngShearVname = "/angle_shear"
MagShearStormVname = "/mag_shear_storm"
AvgMagShearStormVname = "/avg_mag_shear_storm"

# Samples of shear magnitude field at t = 10, 30, 40, 60 hrs sim time.
# These are the beginning and end times of the PSAP and SAP.
MagShear10Vname = "/mag_shear_10"
MagShear20Vname = "/mag_shear_20"
MagShear30Vname = "/mag_shear_30"
MagShear40Vname = "/mag_shear_40"
MagShear50Vname = "/mag_shear_50"
MagShear60Vname = "/mag_shear_60"

# Time step numbering starts with zero
T10 =  20 # time step number at 10 hrs
T20 =  40 # time step number at 20 hrs
T30 =  60 # time step number at 30 hrs
T40 =  80 # time step number at 40 hrs
T50 = 100 # time step number at 50 hrs
T60 = 120 # time step number at 60 hrs

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'  # vertical coords: height
Pname = '/p_coords'  # vertical coords: pressure
Tname = '/t_coords'

# For selection of levels, SFC is layer near surfac 
#              Press coords    Height coords
#  SAL bot        700 mb         3500 m
#  SAL top        500 mb         5500 m
#
#  SFC bot       1000 mb          110 m
#  SFC top        925 mb          761 m

if (Vcoord == 'p'):
    # pressure coords
#    ZsalBot = 700
#    ZsalTop = 500

#    ZsfcBot = 1000
#    ZsfcTop =  925

    ZsalBot = 200
    ZsalTop = 200

    ZsfcBot = 850
    ZsfcTop = 850

    Vunits = "mb"

elif (Vcoord == 'z'):
    # height coords
    ZsalBot = 3000
    ZsalTop = 4000
    
    ZsfcBot =  100
    ZsfcTop =  500

    Vunits = "m"

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Creating wind shear field for simulation: {0:s}".format(Sim))
    print("  Surface layer: Top:    {0:d} {1:s}".format(ZsfcTop, Vunits))
    print("                 Bottom: {0:d} {1:s}".format(ZsfcBot, Vunits))
    print("  SAL:           Top:    {0:d} {1:s}".format(ZsalTop, Vunits))
    print("                 Bottom: {0:d} {1:s}".format(ZsalBot, Vunits))
    print("")

    FilterFname = FilterFnameTemplate.replace("<SIM>", Sim)
    Ufname      = UfnameTemplate.replace("<SIM>", Sim)
    Vfname      = VfnameTemplate.replace("<SIM>", Sim)
    Ofname      = OfnameTemplate.replace("<SIM>", Sim)

    print("  Reading {0:s} ({1:s})".format(FilterFname, FilterVname))
    print("  Reading {0:s} ({1:s})".format(Ufname, Uvname))
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

    # Surface Layer
    if (Vcoord == 'p'):
        Select = np.where((Z <= ZsfcBot) * (Z >= ZsfcTop))
    elif (Vcoord == 'z'):
        Select = np.where((Z >= ZsfcBot) * (Z <= ZsfcTop))

    SFC_Z1 = Select[0][0]
    SFC_Z2 = Select[0][-1]

    # SAL
    if (Vcoord == 'p'):
        Select = np.where((Z <= ZsalBot) * (Z >= ZsalTop))
    elif (Vcoord == 'z'):
        Select = np.where((Z >= ZsalBot) * (Z <= ZsalTop))

    SAL_Z1 = Select[0][0]
    SAL_Z2 = Select[0][-1]

    print(SFC_Z1, SFC_Z2, SAL_Z1, SAL_Z2)

    # Create the output datasets using the COARDS convention so that grads
    # can read these datasets directly (sdfopen).
    UshearDset = h5u.DsetCoards(UshearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    VshearDset = h5u.DsetCoards(VshearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    MagShearDset = h5u.DsetCoards(MagShearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    AngShearDset = h5u.DsetCoards(AngShearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    MagShearStormDset = h5u.DsetCoards(MagShearStormVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    AvgMagShearStormDset = h5u.DsetCoards(AvgMagShearStormVname, 1, [ Nt ])
    MagShear10Dset = h5u.DsetCoards(MagShear10Vname, 2, [ Ny, Nx ])
    MagShear20Dset = h5u.DsetCoards(MagShear20Vname, 2, [ Ny, Nx ])
    MagShear30Dset = h5u.DsetCoards(MagShear30Vname, 2, [ Ny, Nx ])
    MagShear40Dset = h5u.DsetCoards(MagShear40Vname, 2, [ Ny, Nx ])
    MagShear50Dset = h5u.DsetCoards(MagShear50Vname, 2, [ Ny, Nx ])
    MagShear60Dset = h5u.DsetCoards(MagShear60Vname, 2, [ Ny, Nx ])

    UshearDset.Create(Ofile)
    VshearDset.Create(Ofile)
    MagShearDset.Create(Ofile)
    AngShearDset.Create(Ofile)
    MagShearStormDset.Create(Ofile)
    AvgMagShearStormDset.Create(Ofile)
    MagShear10Dset.Create(Ofile)
    MagShear20Dset.Create(Ofile)
    MagShear30Dset.Create(Ofile)
    MagShear40Dset.Create(Ofile)
    MagShear50Dset.Create(Ofile)
    MagShear60Dset.Create(Ofile)

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
        #MAG_SHEAR = np.absolute(U_SHEAR)  # zonal shear
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
   
        # Write out the sampled shear magnitude fields. Add the [...] so
        # that h5py doesn't try to recreate the dataset. We want to
        # create the dataset so that we can attach dimensions. If you
        # simply say "Ofile[MagShear10Vname] = MAG_SHEAR", h5py will try
        # to create the dataset which throws an error since we already
        # created the dataset above.
        if (it == T10):
            Ofile[MagShear10Vname][...] = MAG_SHEAR
        elif (it == T20):
            Ofile[MagShear20Vname][...] = MAG_SHEAR
        elif (it == T30):
            Ofile[MagShear30Vname][...] = MAG_SHEAR
        elif (it == T40):
            Ofile[MagShear40Vname][...] = MAG_SHEAR
        elif (it == T50):
            Ofile[MagShear50Vname][...] = MAG_SHEAR
        elif (it == T60):
            Ofile[MagShear60Vname][...] = MAG_SHEAR

    Ufile.close()
    Vfile.close()

    # Attach dimensions to output fields
    print("  Writing {0:s} ({1:s})".format(Ofname, UshearVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, VshearVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShearVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, AngShearVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShearStormVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, AvgMagShearStormVname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShear10Vname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShear20Vname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShear30Vname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShear40Vname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShear50Vname))
    print("  Writing {0:s} ({1:s})".format(Ofname, MagShear60Vname))

    UshearDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    VshearDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    MagShearDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    AngShearDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    MagShearStormDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    AvgMagShearStormDset.AttachDims(Ofile, Tdim)
    MagShear10Dset.AttachDims(Ofile, Ydim, Xdim)
    MagShear20Dset.AttachDims(Ofile, Ydim, Xdim)
    MagShear30Dset.AttachDims(Ofile, Ydim, Xdim)
    MagShear40Dset.AttachDims(Ofile, Ydim, Xdim)
    MagShear50Dset.AttachDims(Ofile, Ydim, Xdim)
    MagShear60Dset.AttachDims(Ofile, Ydim, Xdim)

    Ofile.close()
    print("")

