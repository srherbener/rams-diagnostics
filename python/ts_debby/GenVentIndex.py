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
import NumUtils as nu

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

WshFnameTemplate = "DIAGS/wind_shear_nv_lite_<SIM>.h5"
WshVname = "/mag_shear"

OutFnameTemplate = "DIAGS/vent_index_<SIM>.h5"
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

WindShearSmoothVname     = "/sm_wind_shear"

SmSatSmoothVname   = "/sm_mid_level_sat_entropy"
SmSmoothVname      = "/sm_mid_level_entropy"
SsstSatSmoothVname = "/sm_sst_sat_entropy"
SbSmoothVname      = "/sm_bl_entropy"
SdefSmoothVname    = "/sm_entropy_deficit"

TsSmoothVname = "/sm_ocean_temp"
ToSmoothVname = "/sm_outflow_temp"
PiSmoothVname = "/sm_pot_intensity"

ViSmoothVname = "/sm_vent_index"

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

    WshFname              = WshFnameTemplate.replace("<SIM>", Sim)
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


    print("  Reading {0:s} ({1:s})".format(WshFname, WshVname))
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
    WshFile              = h5py.File(WshFname, mode='r')
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


    # Create and save the three components of the ventilation index:
    #     wind shear
    #     entropy deficit
    #     potential intensity
    #
    # Use the above components to calculate the ventilation index.
    #
    # Do one time step at a time to save on memory
    #
    WindShear = np.zeros(Nt)

    Sm      = np.zeros(Nt)
    SmSat   = np.zeros(Nt)
    Sb      = np.zeros(Nt)
    SsstSat = np.zeros(Nt)

    Ts = np.zeros(Nt)
    To = np.zeros(Nt)

    for it in np.arange(Nt):
        # read in vars
        MagShear = np.squeeze(WshFile[WshVname][it,...])

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
        Smidlevel    = LayerAverage(S, Dens, Zmid1, Zmid2)
        SsatMidlevel = LayerAverage(Ssat, Dens, Zmid1, Zmid2)
        Sblayer      = LayerAverage(Snv, DensNv, Zbl1, Zbl2)

        Tout = LayerAverage(Tnv, DensNv, Zout1, Zout2)

        ######## Wind Shear ##############
        WindShear[it] = np.mean(MagShear[Radius <= Rout])

        ########## Entropy Deficit ############
        Sm[it]    = np.mean(Smidlevel[(Radius >= Rin) * (Radius <= Rout)])
        SmSat[it] = np.mean(SsatMidlevel[Radius <= Rin])

        Sb[it]      = np.mean(Sblayer[Radius <= Rout])
        SsstSat[it] = np.mean(SsatNvSst[Radius <= Rout])


        ########## Potential Intensity ############
        Ts[it] = np.mean(TsstNv[Radius <= Rout])
        To[it] = np.mean(Tout[Radius <= Rout])


        if ((it % 10) == 0):
            print("  Processing time step: {0:d}".format(it))

    print("")

    # finish calculations
    # entropy deficit
    Sdef = (SmSat - Sm) / (SsstSat - Sb)

    # potiential intensity
    Pi = np.sqrt((Ts - To) * (Ts / To) * CkOverCd * (SsstSat - Sb))

    # ventilation index
    Vi = WindShear * Sdef / Pi

    ##### Smoothed versions ####
    # wind shear
    WindShearSmooth = nu.SmoothLine(WindShear)

    # entropy deficit
    SmSmooth      = nu.SmoothLine(Sm)
    SmSatSmooth   = nu.SmoothLine(SmSat)
    SbSmooth      = nu.SmoothLine(Sb)
    SsstSatSmooth = nu.SmoothLine(SsstSat)

    SdefSmooth = (SmSatSmooth - SmSmooth) / (SsstSatSmooth - Sb)

    # potential intensity
    TsSmooth = nu.SmoothLine(Ts)
    ToSmooth = nu.SmoothLine(To)

    PiSmooth = np.sqrt((TsSmooth - ToSmooth) * (TsSmooth / ToSmooth) * CkOverCd * (SsstSatSmooth - SbSmooth))

    # ventilation index
    ViSmooth = WindShearSmooth * SdefSmooth / PiSmooth


    # Create output with COARDS convention
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

    print("  Writing {0:s} ({1:s})".format(OutFname, WindShearSmoothVname))
    print("")

    print("  Writing {0:s} ({1:s})".format(OutFname, SmSmoothVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SmSatSmoothVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SbSmoothVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SsstSatSmoothVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SdefSmoothVname))
    print("")

    print("  Writing {0:s} ({1:s})".format(OutFname, TsSmoothVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, ToSmoothVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, PiSmoothVname))
    print("")

    print("  Writing {0:s} ({1:s})".format(OutFname, ViSmoothVname))
    print("")

    # Create output with COARDS convention
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

    WindShearSmoothDset = h5u.DsetCoards(WindShearSmoothVname, 1, [ Nt ])

    SmSmoothDset = h5u.DsetCoards(SmSmoothVname, 1, [ Nt ])
    SmSatSmoothDset = h5u.DsetCoards(SmSatSmoothVname, 1, [ Nt ])
    SbSmoothDset = h5u.DsetCoards(SbSmoothVname, 1, [ Nt ])
    SsstSatSmoothDset = h5u.DsetCoards(SsstSatSmoothVname, 1, [ Nt ])
    SdefSmoothDset = h5u.DsetCoards(SdefSmoothVname, 1, [ Nt ])

    TsSmoothDset = h5u.DsetCoards(TsSmoothVname, 1, [ Nt ])
    ToSmoothDset = h5u.DsetCoards(ToSmoothVname, 1, [ Nt ])
    PiSmoothDset = h5u.DsetCoards(PiSmoothVname, 1, [ Nt ])

    ViSmoothDset = h5u.DsetCoards(ViSmoothVname, 1, [ Nt ])


    WindShearDset.Build(OutFile, WindShear, Tdim)

    SmDset.Build(OutFile, Sm, Tdim)
    SmSatDset.Build(OutFile, SmSat, Tdim)
    SbDset.Build(OutFile, Sb, Tdim)
    SsstSatDset.Build(OutFile, SsstSat, Tdim)
    SdefDset.Build(OutFile, Sdef, Tdim)

    TsDset.Build(OutFile, Ts, Tdim)
    ToDset.Build(OutFile, To, Tdim)
    PiDset.Build(OutFile, Pi, Tdim)

    ViDset.Build(OutFile, Vi, Tdim)

    WindShearSmoothDset.Build(OutFile, WindShearSmooth, Tdim)

    SmSmoothDset.Build(OutFile, SmSmooth, Tdim)
    SmSatSmoothDset.Build(OutFile, SmSatSmooth, Tdim)
    SbSmoothDset.Build(OutFile, SbSmooth, Tdim)
    SsstSatSmoothDset.Build(OutFile, SsstSatSmooth, Tdim)
    SdefSmoothDset.Build(OutFile, SdefSmooth, Tdim)

    TsSmoothDset.Build(OutFile, TsSmooth, Tdim)
    ToSmoothDset.Build(OutFile, ToSmooth, Tdim)
    PiSmoothDset.Build(OutFile, PiSmooth, Tdim)

    ViSmoothDset.Build(OutFile, ViSmooth, Tdim)


    # clean up
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
