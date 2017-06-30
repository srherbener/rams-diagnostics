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
StormCenterFnameTemplate = "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
RadiusVname = "/radius"

SstFnameTemplate = "HDF5/<SIM>/HDF5/sst_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SstVname = "/sst"

SpressFnameTemplate = "HDF5/<SIM>/HDF5/sea_press_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SpressVname = "/sea_press"

RmwFnameTemplate = "DIAGS/storm_meas_<SIM>.h5"
RmwVname = "/rmw_t_p"

DensityFnameTemplate = "HDF5/<SIM>/HDF5/density_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
DensityNvFnameTemplate = "HDF5/<SIM>/HDF5/density_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
DensityVname = "/density"

EntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
EntropyNvFnameTemplate = "HDF5/<SIM>/HDF5/entropy_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
EntropyVname = "/entropy"

SatEntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_s_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
SatEntropyNvFnameTemplate = "HDF5/<SIM>/HDF5/entropy_s_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
SatEntropyVname = "/entropy_s"


OutFnameTemplate = "DIAGS/entropy_deficit_<SIM>.h5"
SmSatVname   = "/mid_level_inner_sat_entropy"
SmVname      = "/mid_level_outer_entropy"
SsSatVname   = "/sst_sat_entropy"
SbVname      = "/bl_entropy"
SmlDiffVname = "mid_level_diff_entropy"
SasDiffVname = "sea_air_diff_entropy"
SdVname      = "/entropy_deficit"

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

# For mid-level measurements, use average of 700 to 500 mb levels
# These are indices 3 and 4 along the z-axis
Zmid1 = 3
Zmid2 = 4

# For boundary layer measurements, use average of 1000 to 850 mb levels
# These are indices 0 and 2 along the z-axis
Zbl1 = 0
Zbl2 = 0

# Radius values for averaging regions
Rin  = 150.0 # km
Rout = 300.0 # km

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Calculating entropy deficit for simulation: {0:s}".format(Sim))
    print("")

    StormCenterFname   = StormCenterFnameTemplate.replace("<SIM>", Sim)
    SstFname           = SstFnameTemplate.replace("<SIM>", Sim)
    SpressFname        = SpressFnameTemplate.replace("<SIM>", Sim)
    RmwFname           = RmwFnameTemplate.replace("<SIM>", Sim)
    DensityFname       = DensityFnameTemplate.replace("<SIM>", Sim)
    EntropyFname       = EntropyFnameTemplate.replace("<SIM>", Sim)
    SatEntropyFname    = SatEntropyFnameTemplate.replace("<SIM>", Sim)
    DensityNvFname     = DensityNvFnameTemplate.replace("<SIM>", Sim)
    EntropyNvFname     = EntropyNvFnameTemplate.replace("<SIM>", Sim)
    SatEntropyNvFname  = SatEntropyNvFnameTemplate.replace("<SIM>", Sim)

    OutFname           = OutFnameTemplate.replace("<SIM>", Sim)


    print("  Reading {0:s} ({1:s})".format(StormCenterFname, RadiusVname))
    print("  Reading {0:s} ({1:s})".format(SstFname, SstVname))
    print("  Reading {0:s} ({1:s})".format(SpressFname, SpressVname))
    print("  Reading {0:s} ({1:s})".format(RmwFname, RmwVname))
    print("  Reading {0:s} ({1:s})".format(DensityFname, DensityVname))
    print("  Reading {0:s} ({1:s})".format(EntropyFname, EntropyVname))
    print("  Reading {0:s} ({1:s})".format(SatEntropyFname, SatEntropyVname))
    print("  Reading {0:s} ({1:s})".format(DensityNvFname, DensityVname))
    print("  Reading {0:s} ({1:s})".format(EntropyNvFname, EntropyVname))
    print("  Reading {0:s} ({1:s})".format(SatEntropyNvFname, SatEntropyVname))
    print("")

    # Read coordinates and build the dimensions.
    StormCenterFile   = h5py.File(StormCenterFname, mode='r')
    SstFile           = h5py.File(SstFname, mode='r')
    SpressFile        = h5py.File(SpressFname, mode='r')
    RmwFile           = h5py.File(RmwFname, mode='r')
    DensityFile       = h5py.File(DensityFname, mode='r')
    EntropyFile       = h5py.File(EntropyFname, mode='r')
    SatEntropyFile    = h5py.File(SatEntropyFname, mode='r')
    DensityNvFile     = h5py.File(DensityNvFname, mode='r')
    EntropyNvFile     = h5py.File(EntropyNvFname, mode='r')
    SatEntropyNvFile  = h5py.File(SatEntropyNvFname, mode='r')

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

    SmSat   = np.zeros((Nt))
    Sm      = np.zeros((Nt))
    SsSat   = np.zeros((Nt))
    Sb      = np.zeros((Nt))

    # Entropy deficit is from Tang & Emanuel (2012)
    #
    #     Entropy deficit = (Sm-star - Sm) / (Ssst-star - Sb)
    #
    #     Sm-star: sat. entropy at mid-level taken as average over disk
    #              of radius 150 km, from center of storm
    #     Sm: entropy at mid_level taken as average over annulus
    #              of radius 150 to 300 km, from center of storm
    #     Ssst-star: sat. entropy at SST temperature taken as average
    #              at radius of maximum wind
    #     Sb: entropy in boundary layer taken as average at radius of
    #              maximum wind
    #

    # do one time step at a time to save on memory
    # do formula on columns, then average the results
    for it in np.arange(Nt):
        # read in vars

        # entropy, density at midlevel
        Smidlevel     = np.squeeze(EntropyFile[EntropyVname][it,Zmid1:Zmid2+1,...])
        SsatMidlevel  = np.squeeze(SatEntropyFile[SatEntropyVname][it,Zmid1:Zmid2+1,...])
        DensMidlevel  = np.squeeze(DensityFile[DensityVname][it,Zmid1:Zmid2+1,...])

        # entropy, density at boundary layer
        # use the data with the vortex removed (environmental fields)
        SsatSfc    = np.squeeze(SatEntropyNvFile[SatEntropyVname][it,1,...])
        Sblayer    = np.squeeze(EntropyNvFile[EntropyVname][it,Zbl1:Zbl2+1,...])
        DensBlayer = np.squeeze(DensityNvFile[DensityVname][it,Zbl1:Zbl2+1,...])


        # sea surface temperature (SST) and sea level pressure
        Sst = np.squeeze(SstFile[SstVname][it,...]) + 273.15  # convert to K
        Spress = np.squeeze(SpressFile[SpressVname][it,...]) * 100.0 # covert to Pa

        # radius of maximum wind and radius values
        Rmw = RmwFile[RmwVname][it] / 1000.0 # convert to km
        Radius = np.squeeze(StormCenterFile[RadiusVname][it,...]) # already in km


        # density weighted column averaging
        Smidlevel    = np.squeeze(np.average(Smidlevel, axis=0, weights=DensMidlevel))
        SsatMidlevel = np.squeeze(np.average(SsatMidlevel, axis=0, weights=DensMidlevel))

        #Sblayer = np.squeeze(np.average(Sblayer, axis=0, weights=DensBlayer))


        # Calculate saturation entropy at SST and sea level pressure 
        SsatSst = tu.SatEntropy(Sst, Spress)
        #SsatSst = tu.SatEntropy(Sst, 101300.0)


        # For mid-level, select across the inside radius (disk) for the saturation entropy
        # and across the inside radius to outside radius (annulus) for the outer entropy.
        Sm[it]    = np.mean(Smidlevel[(Radius >= Rin) * (Radius <= Rout)])
        SmSat[it] = np.mean(SsatMidlevel[Radius <= Rin])

        # Select across the entire outside radius disk for both the surface and
        # boundary layer entropies. 
        Sb[it]    = np.mean(Sblayer[(Radius <= 500.0)])
        SsSat[it] = np.mean(SsatSst[Radius <= 500.0])
        #SsSat[it] = np.mean(SsatSfc[Radius <= 500.0])
    

    # Entropy deficit
    SmlDiff = SmSat - Sm
    SasDiff = SsSat - Sb
    Sd = SmlDiff / SasDiff

    # Create output with COARDS convention
    print("  Writing {0:s} ({1:s})".format(OutFname, SmSatVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SmVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SsSatVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SbVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SmlDiffVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SasDiffVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SdVname))

    SmSatDset   = h5u.DsetCoards(SmSatVname, 1, [Nt])
    SmDset      = h5u.DsetCoards(SmVname, 1, [Nt])
    SsSatDset   = h5u.DsetCoards(SsSatVname, 1, [Nt])
    SbDset      = h5u.DsetCoards(SbVname, 1, [Nt])
    SmlDiffDset = h5u.DsetCoards(SmlDiffVname, 1, [Nt])
    SasDiffDset = h5u.DsetCoards(SasDiffVname, 1, [Nt])
    SdDset      = h5u.DsetCoards(SdVname, 1, [Nt])

    SmSatDset.Build(OutFile, SmSat, Tdim)
    SmDset.Build(OutFile, Sm, Tdim)
    SsSatDset.Build(OutFile, SsSat, Tdim)
    SbDset.Build(OutFile, Sb, Tdim)
    SmlDiffDset.Build(OutFile, SmlDiff, Tdim)
    SasDiffDset.Build(OutFile, SasDiff, Tdim)
    SdDset.Build(OutFile, Sd, Tdim)

    print("")

    # clean up
    StormCenterFile.close()
    SstFile.close()
    SpressFile.close()
    RmwFile.close()
    DensityFile.close()
    EntropyFile.close()
    SatEntropyFile.close()
    DensityNvFile.close()
    EntropyNvFile.close()
    SatEntropyNvFile.close()

    OutFile.close()

