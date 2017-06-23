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
Filter150FnameTemplate = "FILTERS/all_150_lite_<SIM>.h5"
Filter150_300FnameTemplate = "FILTERS/all_150_300_lite_<SIM>.h5"
FilterVname = "/filter"

StormCenterFnameTemplate = "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
RadiusVname = "/radius"

SstFnameTemplate = "HDF5/<SIM>/HDF5/sst_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SstVname = "/sst"

SpressFnameTemplate = "HDF5/<SIM>/HDF5/sea_press_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SpressVname = "/sea_press"

RmwFnameTemplate = "DIAGS/storm_meas_<SIM>.h5"
RmwVname = "/rmw_t_p"

DensityFnameTemplate = "HDF5/<SIM>/HDF5/density_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
DensityVname = "/density"

EntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
EntropyVname = "/entropy"

SatEntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_s_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
SatEntropyVname = "/entropy_s"


OutFnameTemplate = "DIAGS/entropy_deficit_<SIM>.h5"
SmSatVname = "/mid_level_inner_sat_entropy"
SmVname    = "/mid_level_outer_entropy"
SsSatVname = "/sst_sat_entropy"
SbVname    = "/bl_entropy"
SdVname    = "/entropy_deficit"

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
Zbl2 = 2

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Calculating entropy deficit for simulation: {0:s}".format(Sim))
    print("")

    Filter150Fname     = Filter150FnameTemplate.replace("<SIM>", Sim)
    Filter150_300Fname = Filter150_300FnameTemplate.replace("<SIM>", Sim)
    StormCenterFname   = StormCenterFnameTemplate.replace("<SIM>", Sim)
    SstFname           = SstFnameTemplate.replace("<SIM>", Sim)
    SpressFname        = SpressFnameTemplate.replace("<SIM>", Sim)
    RmwFname           = RmwFnameTemplate.replace("<SIM>", Sim)
    DensityFname       = DensityFnameTemplate.replace("<SIM>", Sim)
    EntropyFname       = EntropyFnameTemplate.replace("<SIM>", Sim)
    SatEntropyFname    = SatEntropyFnameTemplate.replace("<SIM>", Sim)

    OutFname           = OutFnameTemplate.replace("<SIM>", Sim)


    print("  Reading {0:s} ({1:s})".format(Filter150Fname, FilterVname))
    print("  Reading {0:s} ({1:s})".format(Filter150_300Fname, FilterVname))
    print("  Reading {0:s} ({1:s})".format(StormCenterFname, RadiusVname))
    print("  Reading {0:s} ({1:s})".format(SstFname, SstVname))
    print("  Reading {0:s} ({1:s})".format(SpressFname, SpressVname))
    print("  Reading {0:s} ({1:s})".format(RmwFname, RmwVname))
    print("  Reading {0:s} ({1:s})".format(DensityFname, DensityVname))
    print("  Reading {0:s} ({1:s})".format(EntropyFname, EntropyVname))
    print("  Reading {0:s} ({1:s})".format(SatEntropyFname, SatEntropyVname))
    print("")

    # Read coordinates and build the dimensions.
    Filter150File     = h5py.File(Filter150Fname, mode='r')
    Filter150_300File = h5py.File(Filter150_300Fname, mode='r')
    StormCenterFile   = h5py.File(StormCenterFname, mode='r')
    SstFile           = h5py.File(SstFname, mode='r')
    SpressFile        = h5py.File(SpressFname, mode='r')
    RmwFile           = h5py.File(RmwFname, mode='r')
    DensityFile       = h5py.File(DensityFname, mode='r')
    EntropyFile       = h5py.File(EntropyFname, mode='r')
    SatEntropyFile    = h5py.File(SatEntropyFname, mode='r')

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

    SmSat = np.zeros((Nt))
    Sm    = np.zeros((Nt))
    SsSat = np.zeros((Nt))
    Sb    = np.zeros((Nt))

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
    for it in np.arange(Nt):
        Filter150     = np.squeeze(Filter150File[FilterVname][it,...])
        Filter150_300 = np.squeeze(Filter150_300File[FilterVname][it,...])

        # mid-level -> average 700 and 500mb levels
        S    = np.squeeze(EntropyFile[EntropyVname][it,Zmid1:Zmid2+1,...])
        Ssat = np.squeeze(SatEntropyFile[SatEntropyVname][it,Zmid1:Zmid2+1,...])
        Dens = np.squeeze(DensityFile[DensityVname][it,Zmid1:Zmid2+1,...])

        S    = np.squeeze(np.average(S, axis=0, weights=Dens))
        Ssat = np.squeeze(np.average(Ssat, axis=0, weights=Dens))

        # The filter data has 1's and 0's, so multiply by the corresponding filter and
        # select all non-zero points for the final average value.
        S    = np.multiply(Filter150_300, S)
        Ssat = np.multiply(Filter150, Ssat)

        Sm[it]    = np.mean(S[S > 0.0])
        SmSat[it] = np.mean(Ssat[Ssat > 0.0])

        # boundary layer -> average 1000 to 850 mb levels
        S    = np.squeeze(EntropyFile[EntropyVname][it,Zbl1:Zbl2+1,...])
        Dens = np.squeeze(DensityFile[DensityVname][it,Zbl1:Zbl2+1,...])

        S = np.squeeze(np.average(S, axis=0, weights=Dens))

        # Read the SST and convert this to a saturation entropy based on the pressure
        # being 1013 mb at sea level.
        # storm effects will be ignored since SST is set to a constant field (based on obs)
        Sst = np.squeeze(SstFile[SstVname][it,...])
        Sst = Sst + 273.15
        #Spress = 101300.0
        Spress = np.squeeze(SpressFile[SpressVname][it,...]) * 100.0 # covert to Pa
        Ssat = tu.SatEntropy(Sst, Spress)
    
        # Read in the RMW and radius fields and use these to construct a selection vector for
        # confining the calculations to the region of the RMW.
        Rmw = RmwFile[RmwVname][it] / 1000.0 # convert to km
        Radius = np.squeeze(StormCenterFile[RadiusVname][it,...])
        Radius = np.absolute(Radius - Rmw)
        Select = np.where(Radius <= 20.0)  # select where Radius is within 20 km of RMW

        SsSat[it] = np.mean(Ssat[Select])
        Sb[it] = np.mean(S[Select])

    # Entropy deficit
    Sd = (SmSat - Sm) / (SsSat - Sb)

    # Create output with COARDS convention
    print("  Writing {0:s} ({1:s})".format(OutFname, SmSatVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SmVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SsSatVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SbVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SdVname))

    SmSatDset = h5u.DsetCoards(SmSatVname, 1, [Nt])
    SmDset    = h5u.DsetCoards(SmVname, 1, [Nt])
    SsSatDset = h5u.DsetCoards(SsSatVname, 1, [Nt])
    SbDset    = h5u.DsetCoards(SbVname, 1, [Nt])
    SdDset    = h5u.DsetCoards(SdVname, 1, [Nt])

    SmSatDset.Build(OutFile, SmSat, Tdim)
    SmDset.Build(OutFile, Sm, Tdim)
    SsSatDset.Build(OutFile, SsSat, Tdim)
    SbDset.Build(OutFile, Sb, Tdim)
    SdDset.Build(OutFile, Sd, Tdim)

    print("")

    # clean up
    Filter150File.close()
    Filter150_300File.close()
    StormCenterFile.close()
    SstFile.close()
    SpressFile.close()
    RmwFile.close()
    EntropyFile.close()
    SatEntropyFile.close()

    OutFile.close()

