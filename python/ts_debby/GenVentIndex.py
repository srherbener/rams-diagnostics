#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u
import ThermoUtils as tu

#################################################################################
# SUBROUTINES
#################################################################################

#################################################################################
# LayerAverage()
#
# This routine will do vertical layer averaging. Data is selected from a 3D
# a 3D using two Z values. Vertical averaging is done with density weighting.
#
# Z1, Z2 are index values for the vertical dimension. If they are equal, that
# level is simply selected for the output.
#
def LayerAverage(Var, Dens, Z1, Z2):
    if (Z1 == Z2):
        VarOut = np.squeeze(Var[Z1,:,:])
    else:
        VarOut = np.squeeze(np.average(Var[Z1:Z2+1,:,:], axis = 0, weights=Dens[Z1:Z2+1,:,:]))

    return VarOut


#################################################################################
# MAIN
#################################################################################

Tstring = conf.SetTimeString()

SimList = [
    'TSD_SAL_DUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_NONSAL_NODUST'
    ]
Nsims = len(SimList)

# Use pressure for vertical coordinates (-AP- files). The storm center and SST files come
# in sigma-z coordinates (-AS- files) only, but these are 2D fields meaning that the same
# field would apply to either pressure or sigma-z vertical coordinates.
# 
StormCenterFnameTemplate = "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
RadiusVname = "/radius"

DensityFnameTemplate = "HDF5/<SIM>/HDF5/density_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
DensityNvFnameTemplate = "HDF5/<SIM>/HDF5/density_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
DensityVname = "/density"

EntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
EntropyNvFnameTemplate = "HDF5/<SIM>/HDF5/entropy_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
EntropyVname = "/entropy"

SatEntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_s_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
SatEntropyVname = "/entropy_s"

SstSatEntropyFnameTemplate = "HDF5/<SIM>/HDF5/sst_entropy_s_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SstSatEntropyNvFnameTemplate = "HDF5/<SIM>/HDF5/sst_entropy_s_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SstSatEntropyVname = "/sst_entropy_s"

TempNvFnameTemplate = "HDF5/<SIM>/HDF5/tempc_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
TempVname = "/tempc_basic"

SstNvFnameTemplate = "HDF5/<SIM>/HDF5/sst_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SstVname = "/sst_basic"

UnvFnameTemplate = "HDF5/<SIM>/HDF5/u_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
Uvname = "/u_basic"

VnvFnameTemplate = "HDF5/<SIM>/HDF5/v_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
Vvname = "/v_basic"


OutFnameTemplate = "DIAGS/vent_index_<SIM>.h5"
UshearVname        = "/u_shear"
VshearVname        = "/v_shear"
MagShearVname      = "/mag_shear"
AngShearVname      = "/angle_shear"
MagShearStormVname = "/mag_shear_storm"
WindShearVname     = "/wind_shear"

SmSatVname   = "/mid_level_sat_entropy"
SmVname      = "/mid_level_entropy"
SsstSatVname = "/sst_sat_entropy"
SbVname      = "/bl_entropy"
SdefVname    = "/entropy_deficit"

TsVname = "/ocean_temp"
ToVname = "/outflow_temp"
PiVname = "/pot_intensity"

ViVname = "/vent_index"

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

# For layer averaging, pressure coordinates:
#     level        index
#     1000           0
#      925           1
#      850           2
#      700           3
#      500           4
#      400           5
#      300           6
#      250           7
#      200           8
#      150           9
#      100          10


# For wind shear measurements:
#    surface layer: 1000 - 925 mb
#        SAL layer: 700 - 500 mb
Zsfc1 = 0
Zsfc2 = 1

Zsal1 = 3
Zsal2 = 4

# For mid-level layer: 700 - 500 mb
Zmid1 = 3
Zmid2 = 4

# For boundary layer: 1000 to 850 mb
Zbl1 = 0
Zbl2 = 2

# For outflow layer: 300 - 200 mb
Zout1 = 6
Zout2 = 8

# Radius values for averaging regions
Rin  = 150.0 # km
Rout = 300.0 # km

# Ck/Cd (ratio of enthalpy exchange coefficient to drag coefficient) value is
# from Bell et al. (2012) and Haus et al. (2010). TS Debby was a weak storm
# that ranged from 17 m/s to 22 m/s winds, and in that range the Ck/Cd ratio
# value is roughly 0.6.
CkOverCd = 0.6

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Calculating entropy deficit for simulation: {0:s}".format(Sim))
    print("")

    UnvFname              = UnvFnameTemplate.replace("<SIM>", Sim)
    VnvFname              = VnvFnameTemplate.replace("<SIM>", Sim)
    EntropyFname          = EntropyFnameTemplate.replace("<SIM>", Sim)
    EntropyNvFname        = EntropyNvFnameTemplate.replace("<SIM>", Sim)
    SatEntropyFname       = SatEntropyFnameTemplate.replace("<SIM>", Sim)
    SstSatEntropyNvFname  = SstSatEntropyNvFnameTemplate.replace("<SIM>", Sim)
    DensityFname          = DensityFnameTemplate.replace("<SIM>", Sim)
    DensityNvFname        = DensityNvFnameTemplate.replace("<SIM>", Sim)
    TempNvFname           = TempNvFnameTemplate.replace("<SIM>", Sim)
    SstNvFname            = SstNvFnameTemplate.replace("<SIM>", Sim)
    StormCenterFname      = StormCenterFnameTemplate.replace("<SIM>", Sim)

    OutFname = OutFnameTemplate.replace("<SIM>", Sim)


    print("  Reading {0:s} ({1:s})".format(UnvFname, Uvname))
    print("  Reading {0:s} ({1:s})".format(VnvFname, Uvname))
    print("  Reading {0:s} ({1:s})".format(EntropyFname, EntropyVname))
    print("  Reading {0:s} ({1:s})".format(EntropyNvFname, EntropyVname))
    print("  Reading {0:s} ({1:s})".format(SatEntropyFname, SatEntropyVname))
    print("  Reading {0:s} ({1:s})".format(SstSatEntropyNvFname, SstSatEntropyVname))
    print("  Reading {0:s} ({1:s})".format(DensityFname, DensityVname))
    print("  Reading {0:s} ({1:s})".format(DensityNvFname, DensityVname))
    print("  Reading {0:s} ({1:s})".format(TempNvFname, TempVname))
    print("  Reading {0:s} ({1:s})".format(SstNvFname, SstVname))
    print("  Reading {0:s} ({1:s})".format(StormCenterFname, RadiusVname))
    print("")

    # Read coordinates and build the dimensions.
    UnvFile              = h5py.File(UnvFname, mode='r')
    VnvFile              = h5py.File(VnvFname, mode='r')
    EntropyFile          = h5py.File(EntropyFname, mode='r')
    EntropyNvFile        = h5py.File(EntropyNvFname, mode='r')
    SatEntropyFile       = h5py.File(SatEntropyFname, mode='r')
    SstSatEntropyNvFile  = h5py.File(SstSatEntropyNvFname, mode='r')
    DensityFile          = h5py.File(DensityFname, mode='r')
    DensityNvFile        = h5py.File(DensityNvFname, mode='r')
    TempNvFile           = h5py.File(TempNvFname, mode='r')
    SstNvFile            = h5py.File(SstNvFname, mode='r')
    StormCenterFile      = h5py.File(StormCenterFname, mode='r')

    OutFile = h5py.File(OutFname, mode='w')

    print("    Reading {0:s} ({1:s})".format(EntropyFname, Xname))
    X = EntropyFile[Xname][...]
    Nx = len(X)
    Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
    Xdim.Build(OutFile, X)

    print("    Reading {0:s} ({1:s})".format(EntropyFname, Yname))
    Y = EntropyFile[Yname][...]
    Ny = len(Y)
    Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
    Ydim.Build(OutFile, Y)

    print("    Reading {0:s} ({1:s})".format(EntropyFname, Zname))
    Z = EntropyFile[Zname][...]
    Nz = len(Z)
    Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
    Zdim.Build(OutFile, Z)

    print("    Reading {0:s} ({1:s})".format(EntropyFname, Tname))
    T = EntropyFile[Tname][...]
    Nt = len(T)
    Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
    Tdim.Build(OutFile, T)

    print("")

    # Set up fields for writing output one time step at a time
    UshearDset = h5u.DsetCoards(UshearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    VshearDset = h5u.DsetCoards(VshearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    MagShearDset = h5u.DsetCoards(MagShearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    AngShearDset = h5u.DsetCoards(AngShearVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    MagShearStormDset = h5u.DsetCoards(MagShearStormVname, 3, [ Nt, Ny, Nx ], chunks=( 1, Ny, Nx ))
    WindShearDset = h5u.DsetCoards(WindShearVname, 1, [ Nt ])

    SmDset = h5u.DsetCoards(SmVname, 1, [ Nt ])
    SmSatDset = h5u.DsetCoards(SmSatVname, 1, [ Nt ])
    SbDset = h5u.DsetCoards(SbVname, 1, [ Nt ])
    SsstSatDset = h5u.DsetCoards(SsstSatVname, 1, [ Nt ])
    SdefDset = h5u.DsetCoards(SdefVname, 1, [ Nt ])

    TsDset = h5u.DsetCoards(TsVname, 1, [ Nt ])
    ToDset = h5u.DsetCoards(ToVname, 1, [ Nt ])
    PiDset = h5u.DsetCoards(PiVname, 1, [ Nt ])

    ViDset = h5u.DsetCoards(ViVname, 1, [ Nt ])


    UshearDset.Create(OutFile)
    VshearDset.Create(OutFile)
    MagShearDset.Create(OutFile)
    AngShearDset.Create(OutFile)
    MagShearStormDset.Create(OutFile)
    WindShearDset.Create(OutFile)

    SmDset.Create(OutFile)
    SmSatDset.Create(OutFile)
    SbDset.Create(OutFile)
    SsstSatDset.Create(OutFile)
    SdefDset.Create(OutFile)

    TsDset.Create(OutFile)
    ToDset.Create(OutFile)
    PiDset.Create(OutFile)

    ViDset.Create(OutFile)


    # Create and save the three components of the ventilation index:
    #     wind shear
    #     entropy deficit
    #     potential intensity
    #
    # Use the above components to calculate the ventilation index.
    #
    # Do one time step at a time to save on memory
    #
    for it in np.arange(Nt):
        # read in vars
        Unv = np.squeeze(UnvFile[Uvname][it,...])
        Vnv = np.squeeze(VnvFile[Vvname][it,...])

        S         = np.squeeze(EntropyFile[EntropyVname][it,...])
        Snv       = np.squeeze(EntropyNvFile[EntropyVname][it,...])
        Ssat      = np.squeeze(SatEntropyFile[SatEntropyVname][it,...])
        SsatNvSst = np.squeeze(SstSatEntropyNvFile[SstSatEntropyVname][it,...])

        Tnv    = np.squeeze(TempNvFile[TempVname][it,...]) + 273.15  # convert to Kelvin
        TsstNv = np.squeeze(SstNvFile[SstVname][it,...]) + 273.15

        Dens   = np.squeeze(DensityFile[DensityVname][it,...])
        DensNv = np.squeeze(DensityNvFile[DensityVname][it,...])

        Radius = np.squeeze(StormCenterFile[RadiusVname][it,...])

        # layer averaging on vars (density weighted)
        Usfc = LayerAverage(Unv, DensNv, Zsfc1, Zsfc2)
        Vsfc = LayerAverage(Vnv, DensNv, Zsfc1, Zsfc2)
        Usal = LayerAverage(Unv, DensNv, Zsal1, Zsal2)
        Vsal = LayerAverage(Vnv, DensNv, Zsal1, Zsal2)

        Smidlevel    = LayerAverage(S, Dens, Zmid1, Zmid2)
        SsatMidlevel = LayerAverage(Ssat, Dens, Zmid1, Zmid2)
        Sblayer      = LayerAverage(Snv, DensNv, Zbl1, Zbl2)

        Tout = LayerAverage(Tnv, DensNv, Zout1, Zout2)

        ######## Wind Shear ##############
        Ushear = Usal - Usfc
        Vshear = Vsal - Vsfc

        MagShear = np.sqrt(np.square(Ushear) + np.square(Vshear))
        AngShear = np.arctan2(Vshear, Ushear)

        # This copies MagShear where Radius <= Rout, and places zeros elsewhere
        MagShearStorm = MagShear * (Radius <= Rout)

        WindShear = np.mean(MagShear[Radius <= Rout])


        ########## Entropy Deficit ############
        Sm    = np.mean(Smidlevel[(Radius >= Rin) * (Radius <= Rout)])
        SmSat = np.mean(SsatMidlevel[Radius <= Rin])

        Sb      = np.mean(Sblayer[Radius <= Rout])
        SsstSat = np.mean(SsatNvSst[Radius <= Rout])

        Sdef = (SmSat - Sm) / (SsstSat - Sb)


        ########## Potential Intensity ############
        Ts = np.mean(TsstNv[Radius <= Rout])
        To = np.mean(Tout[Radius <= Rout])

        Pi = np.sqrt((Ts - To) * (Ts / To) * CkOverCd * (SsstSat - Sb))


        ########## Ventilation Index ##############
        Vi = WindShear * Sdef / Pi



        # Write fields into output file
        OutFile[UshearVname][it,...] = Ushear
        OutFile[VshearVname][it,...] = Vshear
        OutFile[MagShearVname][it,...] = MagShear
        OutFile[AngShearVname][it,...] = AngShear
        OutFile[MagShearStormVname][it,...] = MagShearStorm
        OutFile[WindShearVname][it] = WindShear

        OutFile[SmVname][it] = Sm
        OutFile[SmSatVname][it] = SmSat
        OutFile[SbVname][it] = Sb
        OutFile[SsstSatVname][it] = SsstSat
        OutFile[SdefVname][it] = Sdef

        OutFile[TsVname][it] = Ts
        OutFile[ToVname][it] = To
        OutFile[PiVname][it] = Pi

        OutFile[ViVname][it] = Vi


        if ((it % 10) == 0):
            print("  Processing time step: {0:d}".format(it))

    print("")

    # Create output with COARDS convention
    print("  Writing {0:s} ({1:s})".format(OutFname, UshearVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, VshearVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, MagShearVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, AngShearVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, MagShearStormVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, WindShearVname))
    print("")

    print("  Writing {0:s} ({1:s})".format(OutFname, SmVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SmSatVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SbVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SsstSatVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SdefVname))
    print("")

    print("  Writing {0:s} ({1:s})".format(OutFname, TsVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, ToVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, PiVname))
    print("")

    print("  Writing {0:s} ({1:s})".format(OutFname, ViVname))
    print("")

    UshearDset.AttachDims(OutFile, Tdim, Ydim, Xdim)
    VshearDset.AttachDims(OutFile, Tdim, Ydim, Xdim)
    MagShearDset.AttachDims(OutFile, Tdim, Ydim, Xdim)
    AngShearDset.AttachDims(OutFile, Tdim, Ydim, Xdim)
    MagShearStormDset.AttachDims(OutFile, Tdim, Ydim, Xdim)
    WindShearDset.AttachDims(OutFile, Tdim)

    SmDset.AttachDims(OutFile, Tdim)
    SmSatDset.AttachDims(OutFile, Tdim)
    SbDset.AttachDims(OutFile, Tdim)
    SsstSatDset.AttachDims(OutFile, Tdim)
    SdefDset.AttachDims(OutFile, Tdim)

    TsDset.AttachDims(OutFile, Tdim)
    ToDset.AttachDims(OutFile, Tdim)
    PiDset.AttachDims(OutFile, Tdim)

    ViDset.AttachDims(OutFile, Tdim)


    # clean up
    UnvFile.close()
    VnvFile.close()
    EntropyFile.close()
    EntropyNvFile.close()
    SatEntropyFile.close()
    SstSatEntropyNvFile.close()
    DensityFile.close()
    DensityNvFile.close()
    TempNvFile.close()
    SstNvFile.close()
    StormCenterFile.close()

    OutFile.close()





###         # entropy, density at midlevel
###         Smidlevel     = np.squeeze(EntropyFile[EntropyVname][it,Zmid1:Zmid2+1,...])
###         SsatMidlevel  = np.squeeze(SatEntropyFile[SatEntropyVname][it,Zmid1:Zmid2+1,...])
###         DensMidlevel  = np.squeeze(DensityFile[DensityVname][it,Zmid1:Zmid2+1,...])
### 
###         # entropy, density at boundary layer
###         # use the data with the vortex removed (environmental fields)
###         SsatSfc    = np.squeeze(SatEntropyNvFile[SatEntropyVname][it,1,...])
###         Sblayer    = np.squeeze(EntropyNvFile[EntropyVname][it,Zbl1:Zbl2+1,...])
###         DensBlayer = np.squeeze(DensityNvFile[DensityVname][it,Zbl1:Zbl2+1,...])
### 
###         # saturation entropy at SST values
###         SsatSst = np.squeeze(SstSatEntropyNvFile[SstSatEntropyVname][it,...])
### 
###         # radius of maximum wind and radius values
###         Radius = np.squeeze(StormCenterFile[RadiusVname][it,...]) # already in km
### 
### 
###         # density weighted column averaging
###         Smidlevel    = np.squeeze(np.average(Smidlevel, axis=0, weights=DensMidlevel))
###         SsatMidlevel = np.squeeze(np.average(SsatMidlevel, axis=0, weights=DensMidlevel))
### 
###         Sblayer = np.squeeze(np.average(Sblayer, axis=0, weights=DensBlayer))
### 
### 
###         # For mid-level, select across the inside radius (disk) for the saturation entropy
###         # and across the inside radius to outside radius (annulus) for the outer entropy.
###         Sm[it]    = np.mean(Smidlevel[(Radius >= Rin) * (Radius <= Rout)])
###         SmSat[it] = np.mean(SsatMidlevel[Radius <= Rin])
### 
###         # For surface, select across the entire storm location from center to outside radius.
###         SbSel     = Sblayer[(Radius <= Rout)]  # yeilds a linear array with selected points
###         SstSatSel = SsatSst[(Radius <= Rout)]
### 
###         Sb[it]    = np.mean(SbSel)
###         SsstSat[it] = np.mean(SstSatSel)
### 
###         #SasDiff[it] = np.mean(SstSatSel - SbSel)  # average of differences
### 
###     # Entropy deficit
###     SmlDiff = SmSat - Sm
###     SasDiff = SsstSat - Sb  # difference of averages
###     Sd = SmlDiff / SasDiff
### 
