#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))

import h5py
import numpy as np
import ConfigRce as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()

SimList = [
    'RCE_3km_1mom',
    'RCE_3km_1mom_db',
    'RCE_3km_1mom_db_udef',
    'RCE_3km_1mom_db_rlongup',
    'RCE_3km_1mom_dm',
#    'RCE_3km_2mom',
    'RCE_3km_2mom_db',
    'RCE_3km_2mom_db_udef',
    'RCE_3km_2mom_db_rlongup',
    'RCE_3km_2mom_dm',
    'RCE_3km_2mom_dm_lrz',
    ]
Nsims = len(SimList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

InFileTemplate = "HDF5/<SIM>/HDF5/<FPREFIX>-<SIM>-LS-2012-01-01-000000-g1.h5"

def WriteOutputDset(OutFname, OutVname, Nt, OutVar, Tdim):
    print("  Writing: {0:s} ({1:s})".format(OutFname, OutVname))
    Dset = h5u.DsetCoards(OutVname, 1, [ Nt ])
    Dset.Build(Ofile, OutVar, Tdim)

for isim in range(Nsims):
    Sim = SimList[isim]
    print("************************************************************************************")
    print("Generating equilibrium measurements for simulation: {0:s}".format(Sim))
    print('')

    # Form input file names
    LatFluxFname  = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "lat_flux")
    SensFluxFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "sens_flux")

    SfcSwdnFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "swdn_sfc")
    SfcAlbFname  = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "albedt")
    SfcLwdnFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "lwdn_sfc")
    SfcLwupFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "lwup_sfc")

    TopSwdnFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "swdn_toa")
    TopSwupFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "swup_toa")
    TopLwupFname = InFileTemplate.replace("<SIM>", Sim).replace("<FPREFIX>", "lwup_toa")

    # Corresponding input file dataset names
    LatFluxVname = "/lat_flux"
    SensFluxVname = "/sens_flux"

    SfcSwdnVname = "/rshort"
    SfcAlbVname  = "/albedt"
    SfcLwdnVname = "/rlong"
    SfcLwupVname = "/rlongup"

    TopSwdnVname = "/swdn"
    TopSwupVname = "/swup"
    TopLwupVname = "/lwup"

    # Open the input files
    LatFluxFile  = h5py.File(LatFluxFname, mode='r')
    SensFluxFile = h5py.File(SensFluxFname, mode='r')

    SfcSwdnFile = h5py.File(SfcSwdnFname, mode='r')
    SfcAlbFile = h5py.File(SfcAlbFname, mode='r')
    SfcLwdnFile = h5py.File(SfcLwdnFname, mode='r')
    SfcLwupFile = h5py.File(SfcLwupFname, mode='r')

    TopSwdnFile = h5py.File(TopSwdnFname, mode='r')
    TopSwupFile = h5py.File(TopSwupFname, mode='r')
    TopLwupFile = h5py.File(TopLwupFname, mode='r')

    # Read in the coordinates from the latent heat flux file
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, Xname))
    X = LatFluxFile[Xname][...]
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, Yname))
    Y = LatFluxFile[Yname][...]
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, Zname))
    Z = LatFluxFile[Zname][...]
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, Tname))
    T = LatFluxFile[Tname][...]
    print("")

    Nx = len(X)
    Ny = len(Y)
    Nz = len(Z)
    Nt = len(T)

    # Inform user that we are reading the input files
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, LatFluxVname))
    print("  Reading: {0:s} ({1:s})".format(SensFluxFname, SensFluxVname))

    print("  Reading: {0:s} ({1:s})".format(SfcSwdnFname, SfcSwdnVname))
    print("  Reading: {0:s} ({1:s})".format(SfcAlbFname, SfcAlbVname))
    print("  Reading: {0:s} ({1:s})".format(SfcLwdnFname, SfcLwdnVname))
    print("  Reading: {0:s} ({1:s})".format(SfcLwupFname, SfcLwupVname))

    print("  Reading: {0:s} ({1:s})".format(TopSwdnFname, TopSwdnVname))
    print("  Reading: {0:s} ({1:s})".format(TopSwupFname, TopSwupVname))
    print("  Reading: {0:s} ({1:s})".format(TopLwupFname, TopLwupVname))

    print("")

    # Open the output file
    OutFname = "DIAGS/eq_meas_{0:s}.h5".format(Sim)
    Ofile = h5py.File(OutFname, mode='w')

    # Write out the coordinate data, mark these as dimensions
    Xdim = h5u.DimCoards(Xname, 1, [ Nx ], 'x')
    Xdim.Build(Ofile, X)

    Ydim = h5u.DimCoards(Yname, 1, [ Ny ], 'y')
    Ydim.Build(Ofile, Y)

    Zdim = h5u.DimCoards(Zname, 1, [ Nz ], 'z')
    Zdim.Build(Ofile, Z)

    Tdim = h5u.DimCoards(Tname, 1, [ Nt ], 't', tstring=Tstring)
    Tdim.Build(Ofile, T)

    # Process one time step at a time
    AvgSfcLat  = np.zeros(Nt)
    AvgSfcSens = np.zeros(Nt)

    AvgSfcSwdn = np.zeros(Nt)
    AvgSfcSwup = np.zeros(Nt)
    AvgSfcLwdn = np.zeros(Nt)
    AvgSfcLwup = np.zeros(Nt)

    AvgTopSwdn = np.zeros(Nt)
    AvgTopSwup = np.zeros(Nt)
    AvgTopLwup = np.zeros(Nt)

    Thf  = np.zeros(Nt)
    Qrad = np.zeros(Nt)

    for i in range(Nt):
        # Exclude the horizontal boundaries since RAMS uses those for boundary
        # conditions.
        SfcLat = np.squeeze(LatFluxFile[LatFluxVname][i,...])[1:-1,1:-1]
        SfcSens = np.squeeze(SensFluxFile[SensFluxVname][i,...])[1:-1,1:-1]

        SfcSwdn = np.squeeze(SfcSwdnFile[SfcSwdnVname][i,...])[1:-1,1:-1]
        SfcAlb  = np.squeeze(SfcAlbFile[SfcAlbVname][i,...])[1:-1,1:-1]
        SfcLwdn = np.squeeze(SfcLwdnFile[SfcLwdnVname][i,...])[1:-1,1:-1]
        SfcLwup = np.squeeze(SfcLwupFile[SfcLwupVname][i,...])[1:-1,1:-1]

        TopSwdn = np.squeeze(TopSwdnFile[TopSwdnVname][i,...])[1:-1,1:-1]
        TopSwup = np.squeeze(TopSwupFile[TopSwupVname][i,...])[1:-1,1:-1]
        TopLwup = np.squeeze(TopLwupFile[TopLwupVname][i,...])[1:-1,1:-1]

        SfcSwup = SfcSwdn * SfcAlb

        # Form Thf and Qrad
        #  Thf is the sum of latent heat and sensible heat fluxes
        #
        #  Qrad is the net radiative flux of the entire column
        #    or: (Net upward flux at TOA) - (Net upward flux at sfc)
        #
        #    where:
        #      (Net upward flux at TOA) = (TopSwup + TopLwup) - TopSwdn
        #      (Net upward flux at sfc) = (SfcSwup + SfcLwup) - (SfcSwdn + SfcLwdn)
        #
        #  Calculate Thf and Qrad per column, then take average
        #
        #  Save out averages of all teh component quantities too.
        #
        # Use double precision since the size of the arrays can get large.
        TempVar = SfcLat + SfcSens
        Thf[i] = np.mean(TempVar.astype(np.float64))

        TempVar = ((TopSwup + TopLwup) - (TopSwdn)) - ((SfcSwup + SfcLwup) - (SfcSwdn + SfcLwdn))
        Qrad[i] = np.mean(TempVar.astype(np.float64))

        AvgSfcLat[i] = np.mean(SfcLat)
        AvgSfcSens[i] = np.mean(SfcSens)

        AvgSfcSwdn[i] = np.mean(SfcSwdn)
        AvgSfcSwup[i] = np.mean(SfcSwup)
        AvgSfcLwdn[i] = np.mean(SfcLwdn)
        AvgSfcLwup[i] = np.mean(SfcLwup)

        AvgTopSwdn[i] = np.mean(TopSwdn)
        AvgTopSwup[i] = np.mean(TopSwup)
        AvgTopLwup[i] = np.mean(TopLwup)

        if ((i % 10) == 0):
            print("  Processing time step: {0:d}".format(i))

    print("")

    # close input files
    LatFluxFile.close()
    SensFluxFile.close()

    SfcSwdnFile.close()
    SfcAlbFile.close()
    SfcLwdnFile.close()
    SfcLwupFile.close()

    TopSwdnFile.close()
    TopSwupFile.close()
    TopLwupFile.close()

    # Build the output datasets
    WriteOutputDset(OutFname, "/therm_heat_flux", Nt, Thf, Tdim);
    WriteOutputDset(OutFname, "/rad_flux_div", Nt, Qrad, Tdim);

    WriteOutputDset(OutFname, "/sfc_lat", Nt, AvgSfcLat, Tdim);
    WriteOutputDset(OutFname, "/sfc_sens", Nt, AvgSfcSens, Tdim);

    WriteOutputDset(OutFname, "/sfc_swdn", Nt, AvgSfcSwdn, Tdim);
    WriteOutputDset(OutFname, "/sfc_swup", Nt, AvgSfcSwup, Tdim);
    WriteOutputDset(OutFname, "/sfc_lwdn", Nt, AvgSfcLwdn, Tdim);
    WriteOutputDset(OutFname, "/sfc_lwup", Nt, AvgSfcLwup, Tdim);

    WriteOutputDset(OutFname, "/top_swdn", Nt, AvgTopSwdn, Tdim);
    WriteOutputDset(OutFname, "/top_swup", Nt, AvgTopSwup, Tdim);
    WriteOutputDset(OutFname, "/top_lwup", Nt, AvgTopLwup, Tdim);

    print("")

    Ofile.close()
