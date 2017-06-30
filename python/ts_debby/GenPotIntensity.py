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
# Use the vortex subtracted quantities since in Tang and Emanuel (2012), the authors wanted
# to remove the effects of the storm on the potential intensity calculation.
StormCenterFnameTemplate = "HDF5/<SIM>/HDF5/storm_center_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
RadiusVname = "/radius"

TempFnameTemplate = "HDF5/<SIM>/HDF5/tempc_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
TempVname = "/tempc"

SstFnameTemplate = "HDF5/<SIM>/HDF5/sst_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SstVname = "/sst"

SpressFnameTemplate = "HDF5/<SIM>/HDF5/sea_press_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SpressVname = "/sea_press"

RmwFnameTemplate = "DIAGS/storm_meas_<SIM>.h5"
RmwVname = "/rmw_t_p"

DensityFnameTemplate = "HDF5/<SIM>/HDF5/density_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
DensityVname = "/density"

EntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
EntropyVname = "/entropy"

SatEntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_s_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
SatEntropyVname = "/entropy_s"


OutFnameTemplate = "DIAGS/pot_intensity_<SIM>.h5"
SsSatVname = "/sst_sat_entropy"
SbVname    = "/bl_entropy"
PiVname    = "/pot_intensity"

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

# Ck/Cd (ratio of enthalpy exchange coefficient to drag coefficient) value is
# from Bell et al. (2012) and Haus et al. (2010). TS Debby was a weak storm
# that ranged from 17 m/s to 22 m/s winds, and in that range the Ck/Cd ratio
# value is roughly 0.6.
CkOverCd = 0.6
 
# For boundary layer measurements, use average of 1000 to 850 mb levels
# These are indices 0 and 2 along the z-axis
Zbl1 = 0
Zbl2 = 2

# For the outflow level (~9km) use 300mb which is index 6
Zout = 6

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Calculating entropy deficit for simulation: {0:s}".format(Sim))
    print("")

    StormCenterFname   = StormCenterFnameTemplate.replace("<SIM>", Sim)
    TempFname          = TempFnameTemplate.replace("<SIM>", Sim)
    SstFname           = SstFnameTemplate.replace("<SIM>", Sim)
    SpressFname        = SpressFnameTemplate.replace("<SIM>", Sim)
    RmwFname           = RmwFnameTemplate.replace("<SIM>", Sim)
    DensityFname       = DensityFnameTemplate.replace("<SIM>", Sim)
    EntropyFname       = EntropyFnameTemplate.replace("<SIM>", Sim)
    SatEntropyFname    = SatEntropyFnameTemplate.replace("<SIM>", Sim)

    OutFname           = OutFnameTemplate.replace("<SIM>", Sim)


    print("  Reading {0:s} ({1:s})".format(StormCenterFname, RadiusVname))
    print("  Reading {0:s} ({1:s})".format(TempFname, TempVname))
    print("  Reading {0:s} ({1:s})".format(SstFname, SstVname))
    print("  Reading {0:s} ({1:s})".format(SpressFname, SpressVname))
    print("  Reading {0:s} ({1:s})".format(RmwFname, RmwVname))
    print("  Reading {0:s} ({1:s})".format(DensityFname, DensityVname))
    print("  Reading {0:s} ({1:s})".format(EntropyFname, EntropyVname))
    print("  Reading {0:s} ({1:s})".format(SatEntropyFname, SatEntropyVname))
    print("")

    # Read coordinates and build the dimensions.
    StormCenterFile   = h5py.File(StormCenterFname, mode='r')
    TempFile          = h5py.File(TempFname, mode='r')
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

    SsSat = np.zeros((Nt))
    Sb    = np.zeros((Nt))

    # From Tang & Emanuel (2012)
    #
    #    Vm^2 = (Ts *(Ts-To_bar))/To_bar * CkOverCd * (Ssst-star -  Sb)
    #
    #    To is from outflow region. In these simulations, the outflow is
    #       at 9km elevation or roughly 300 mb. To_bar is mean outflow
    #       temperature.
    #    Ts is sst
    #
    #    air-sea disequilibrium = (Ssst-star - Sb)
    #       Ssst-star: sat. entropy at SST temperature taken as average
    #              at radius of maximum wind
    #       Sb: entropy in boundary layer taken as average at radius of
    #              maximum wind
    #

    # do one time step at a time to save on memory
    for it in np.arange(Nt):
        Tout = np.squeeze(TempFile[TempVname][it,Zout,...]) + 273.15 # convert to K

        # boundary layer -> average 1000 to 850 mb levels
        S    = np.squeeze(EntropyFile[EntropyVname][it,Zbl1:Zbl2+1,...])
        Dens = np.squeeze(DensityFile[DensityVname][it,Zbl1:Zbl2+1,...])

        S = np.squeeze(np.average(S, axis=0, weights=Dens))

        # Read the SST and convert this to a saturation entropy based on the pressure
        # being 1013 mb at sea level.
        # storm effects will be ignored since SST is set to a constant field (based on obs)
        Sst = np.squeeze(SstFile[SstVname][it,...])
        Sst = Sst + 273.15  # convert to K
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

    # Air-sea disequilibrium
    ASdeq = SsSat - Sb

    # Create output with COARDS convention
    print("  Writing {0:s} ({1:s})".format(OutFname, SsSatVname))
    print("  Writing {0:s} ({1:s})".format(OutFname, SbVname))

    SsSatDset = h5u.DsetCoards(SsSatVname, 1, [Nt])
    SbDset    = h5u.DsetCoards(SbVname, 1, [Nt])

    SsSatDset.Build(OutFile, SsSat, Tdim)
    SbDset.Build(OutFile, Sb, Tdim)

    print("")

    # clean up
    StormCenterFile.close()
    SstFile.close()
    SpressFile.close()
    RmwFile.close()
    EntropyFile.close()
    SatEntropyFile.close()

    OutFile.close()

