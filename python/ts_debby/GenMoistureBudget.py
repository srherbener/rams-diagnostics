#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))
import re

import h5py
import numpy as np
import ConfigTsd as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()
RhoAir = conf.SetRhoAir()

SimList = [
    'TSD_NONSAL_NODUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_SAL_DUST',
    ]
Nsims = len(SimList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

# Input file paths, dataset names
LatFluxFtemp = "DIAGS/hda_meas_ts_lat_flux_<SIM>.h5"
#LatFluxVname = "/nw_lat_flux"
LatFluxVname = "/nw_nstorm_lat_flux"
VintVaporFtemp  = "DIAGS/hda_meas_ts_vint_vapor_<SIM>.h5"
#VintVaporVname  = "/nw_vint_vapor"
VintVaporVname  = "/nw_nstorm_vint_vapor"
TotPrecipFtemp = "DIAGS/hda_meas_ts_precip_<SIM>.h5"
#TotPrecipVname = "/nw_accpcp_mass"
TotPrecipVname = "/nw_nstorm_accpcp_mass"

# Output dataset names
SfcMoistVname = "/sfc_moisture"
PwaterVname = "/pwater"
PrecipVname = "/precip"
AdvectVname = "adv_in_moisture"

Lv = 2.5e6  # J/kg, latent heat of vaporization for water

TimeInterval = 1800.0 # sec, 1/2 hour between time steps

for isim in range(Nsims):
    Sim = SimList[isim]

    print("************************************************************************************")
    print("Generating rain drop energy for simulation: {0:s}".format(Sim))
    print('')

    # Form input file names
    LatFluxFname = LatFluxFtemp.replace("<SIM>", Sim)
    VintVaporFname = VintVaporFtemp.replace("<SIM>",  Sim)
    TotPrecipFname = TotPrecipFtemp.replace("<SIM>", Sim)

    # Open the input files
    LatFluxFile = h5py.File(LatFluxFname, mode='r')
    VintVaporFile = h5py.File(VintVaporFname,  mode='r')
    TotPrecipFile = h5py.File(TotPrecipFname, mode='r')

    # Read in the coordinates from the temperature difference file
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, Xname))
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, Yname))
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, Zname))
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, Tname))
    print("")

    X = LatFluxFile[Xname][...]
    Y = LatFluxFile[Yname][...]
    Z = LatFluxFile[Zname][...]
    T = LatFluxFile[Tname][...]

    Nx = len(X)
    Ny = len(Y)
    Nz = len(Z)
    Nt = len(T)

    # Inform user that we are reading the input files
    print("  Reading: {0:s} ({1:s})".format(LatFluxFname, LatFluxVname))
    print("  Reading: {0:s} ({1:s})".format(VintVaporFname, VintVaporVname))
    print("  Reading: {0:s} ({1:s})".format(TotPrecipFname, TotPrecipVname))
    print("")

    # Open the output file
    OutFname = "DIAGS/moisture_budget_{0:s}.h5".format(Sim)
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

    # Build the output datasets:
    SfcMoistDset = h5u.DsetCoards(SfcMoistVname, 1, [Nt], chunks=())
    SfcMoistDset.Create(Ofile)

    PwaterDset = h5u.DsetCoards(PwaterVname, 1, [Nt], chunks=())
    PwaterDset.Create(Ofile)

    PrecipDset = h5u.DsetCoards(PrecipVname, 1, [Nt], chunks=())
    PrecipDset.Create(Ofile)

    AdvectDset = h5u.DsetCoards(AdvectVname, 1, [Nt], chunks=())
    AdvectDset.Create(Ofile)

    # read in the variables and construct the moisture budget
    LatFlux = np.squeeze(LatFluxFile[LatFluxVname][...])        # W/m^2
    VintVapor = np.squeeze(VintVaporFile[VintVaporVname][...])  # kg/m^2
    TotPrecip = np.squeeze(TotPrecipFile[TotPrecipVname][...])  # kg/m^2

    # close input files
    LatFluxFile.close()
    VintVaporFile.close()
    TotPrecipFile.close()

    # The moisture budget is given by:
    #    PW = AdvectMoistIn + SfcMoist - Precip
    #
    #    AdvectMoistIn can be positive (moisture coming into the region)
    #    or negative (moisture leaving the region)

    # These go into the budget as is
    Pwater = VintVapor
    Precip = TotPrecip

    # The latent heat flux can be converted to a moisture flux by dividing
    # by Lv (latent heat of vaporization, J/kg). This yields kg/m^2/s which
    # is an instantaneous rate (surface moisture flux). There is no way of
    # knowing what happened in between time steps, yet we need to
    # reconstruct the accumulation of moisture transferred from the surface
    # into the atmosphere.
    #
    # Make the assumption that the surface moisture flux (rate) is smoothly
    # changing, and the 30 minute interval between time steps is short
    # enough that we are not missing any periodic fluctuations (such as
    # a diurnal cycle).
    #
    # Then form a representative rate between two time steps that from
    # the average rate of the instantaneous rates at the two time steps.
    # This representative rate can be multiplied by the time interval to
    # yeild the incremental amount of moisture that was transferred into
    # the atmosphere. Form the starting point, from the initial values of
    # the precipitable water and precip amounts (assume that the contribution
    # of advection is zero at the initial point).
    SfcMoistFlux = (LatFlux / Lv )   # kg/m^2/s
    AvgSfcMoistFlux = (SfcMoistFlux[:-1] + SfcMoistFlux[1:]) * 0.5

    SfcMoist = np.zeros([Nt])
    SfcMoist[0] = Pwater[0] + Precip[0]
    for i in range(1,Nt):
        # AvgSfcMoistFlux is length Nt-1, where the n-th element is the
        # representative rate between SfcMoist[n-1] and SfcMoist[n]. 
        SfcMoist[i] = SfcMoist[i-1] + (AvgSfcMoistFlux[i-1] * TimeInterval)

    # The advective piece is formed as a residual in the budget.
    AdvectMoistIn = Pwater - SfcMoist + Precip


    Ofile[SfcMoistVname][...] = SfcMoist
    Ofile[PwaterVname][...] = Pwater
    Ofile[PrecipVname][...] = Precip
    Ofile[AdvectVname][...] = AdvectMoistIn

    # Attach dimensions to output datasets
    print("  Writing {0:s} ({1:s})".format(OutFname, SfcMoistVname))
    SfcMoistDset.AttachDims(Ofile, Tdim)
    print("  Writing {0:s} ({1:s})".format(OutFname, PwaterVname))
    PwaterDset.AttachDims(Ofile, Tdim)
    print("  Writing {0:s} ({1:s})".format(OutFname, PrecipVname))
    PrecipDset.AttachDims(Ofile, Tdim)
    print("  Writing {0:s} ({1:s})".format(OutFname, AdvectVname))
    AdvectDset.AttachDims(Ofile, Tdim)

    print("")

    Ofile.close()
