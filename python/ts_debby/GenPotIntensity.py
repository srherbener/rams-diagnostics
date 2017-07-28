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

TempFnameTemplate = "HDF5/<SIM>/HDF5/tempc_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
TempVname = "/tempc_basic"

SstFnameTemplate = "HDF5/<SIM>/HDF5/sst_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SstVname = "/sst_basic"

DensityFnameTemplate = "HDF5/<SIM>/HDF5/density_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
DensityVname = "/density"

EntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
EntropyVname = "/entropy"

SatEntropyFnameTemplate = "HDF5/<SIM>/HDF5/entropy_s_nv_lite-<SIM>-AP-2006-08-20-120000-g3.h5"
SatEntropyVname = "/entropy_s"

SstSatEntropyFnameTemplate = "HDF5/<SIM>/HDF5/sst_entropy_s_nv_lite-<SIM>-AS-2006-08-20-120000-g3.h5"
SstSatEntropyVname = "/sst_entropy_s"


OutFnameTemplate = "DIAGS/pot_intensity_<SIM>.h5"
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

# Radius values for the averaging
Rout = 300.0 # km

for isim in range(Nsims):
    Sim = SimList[isim]
    print("*****************************************************************")
    print("Calculating potential intensity for simulation: {0:s}".format(Sim))
    print("")

    StormCenterFname   = StormCenterFnameTemplate.replace("<SIM>", Sim)
    TempFname          = TempFnameTemplate.replace("<SIM>", Sim)
    SstFname           = SstFnameTemplate.replace("<SIM>", Sim)
    DensityFname       = DensityFnameTemplate.replace("<SIM>", Sim)
    EntropyFname       = EntropyFnameTemplate.replace("<SIM>", Sim)
    SatEntropyFname    = SatEntropyFnameTemplate.replace("<SIM>", Sim)
    SstSatEntropyFname = SstSatEntropyFnameTemplate.replace("<SIM>", Sim)

    OutFname           = OutFnameTemplate.replace("<SIM>", Sim)


    print("  Reading {0:s} ({1:s})".format(StormCenterFname, RadiusVname))
    print("  Reading {0:s} ({1:s})".format(TempFname, TempVname))
    print("  Reading {0:s} ({1:s})".format(SstFname, SstVname))
    print("  Reading {0:s} ({1:s})".format(DensityFname, DensityVname))
    print("  Reading {0:s} ({1:s})".format(EntropyFname, EntropyVname))
    print("  Reading {0:s} ({1:s})".format(SatEntropyFname, SatEntropyVname))
    print("  Reading {0:s} ({1:s})".format(SstSatEntropyFname, SstSatEntropyVname))
    print("")

    # Read coordinates and build the dimensions.
    StormCenterFile   = h5py.File(StormCenterFname, mode='r')
    TempFile          = h5py.File(TempFname, mode='r')
    SstFile           = h5py.File(SstFname, mode='r')
    DensityFile       = h5py.File(DensityFname, mode='r')
    EntropyFile       = h5py.File(EntropyFname, mode='r')
    SatEntropyFile    = h5py.File(SatEntropyFname, mode='r')
    SstSatEntropyFile = h5py.File(SstSatEntropyFname, mode='r')

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

    SsSat   = np.zeros((Nt))
    Sb      = np.zeros((Nt))
    Ts      = np.zeros((Nt))
    To      = np.zeros((Nt))
    SasDiff = np.zeros((Nt))

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
        # ocean temperature
        Tsst = np.squeeze(SstFile[SstVname][it,...]) + 273.15 # convert to K

        # temperature at outflow
        Tout = np.squeeze(TempFile[TempVname][it,Zout,...]) + 273.15 # convert to K

        # entropy, density at boundary layer
        # use the data with the vortex removed (environmental fields)
        Sblayer    = np.squeeze(EntropyFile[EntropyVname][it,Zbl1:Zbl2+1,...])
        DensBlayer = np.squeeze(DensityFile[DensityVname][it,Zbl1:Zbl2+1,...])

        # saturation entropy at SST values
        SsatSst = np.squeeze(SstSatEntropyFile[SstSatEntropyVname][it,...])

        # radius of maximum wind and radius values
        Radius = np.squeeze(StormCenterFile[RadiusVname][it,...]) # already in km

        # density weighted column averaging
        Sblayer = np.squeeze(np.average(Sblayer, axis=0, weights=DensBlayer))

        # For surface, select across the entire storm location from center to outside radius.
        SbSel     = Sblayer[(Radius <= Rout)]  # yeilds a linear array with selected points
        SstSatSel = SsatSst[(Radius <= Rout)]

        Sb[it]    = np.mean(SbSel)
        SsSat[it] = np.mean(SstSatSel)

        #SasDiff[it] = np.mean(SstSatSel - SbSel)  # average of differences

        Ts[it] = np.mean(Tsst[(Radius <= Rout)])
        To[it] = np.mean(Tout[(Radius <= Rout)])

    # Potential Intensity
    SasDiff = SsSat - Sb # difference of averages

    PotInt =  np.sqrt((Ts - To) * (Ts / To) * CkOverCd * SasDiff)

    # Create output with COARDS convention
    print("  Writing {0:s} ({1:s})".format(OutFname, PiVname))

    PiDset = h5u.DsetCoards(PiVname, 1, [Nt])

    PiDset.Build(OutFile, PotInt, Tdim)

    print("")

    # clean up
    StormCenterFile.close()
    TempFile.close()
    DensityFile.close()
    EntropyFile.close()
    SatEntropyFile.close()
    SstSatEntropyFile.close()

    OutFile.close()

