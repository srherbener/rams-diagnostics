#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))
import re

import h5py
import numpy as np
import ConfigRce as conf
import Hdf5Utils as h5u


Tstring = conf.SetTimeString()

SimList = [
    'RCE_3km_1mom',
    'RCE_3km_1mom_db',
    'RCE_3km_1mom_dm',
#    'RCE_3km_2mom',
    'RCE_3km_2mom_db',
    'RCE_3km_2mom_dm',
    'RCE_3km_2mom_dm_lrz',
    ]
Nsims = len(SimList)

Xname = '/x_coords'
Yname = '/y_coords'
Zname = '/z_coords'
Tname = '/t_coords'

# Input file paths, dataset names
TdiffFtemp = "HDF5/<SIM>/HDF5/rain_air_tempdif_sfc-<SIM>-AS-2012-01-01-000000-g1.h5"
TdiffVname = "/rain_air_tempdif"
DensFtemp  = "HDF5/<SIM>/HDF5/dn0_sfc-<SIM>-AS-2012-01-01-000000-g1.h5"
DensVname  = "/dn0"
RmassFtemp = "HDF5/<SIM>/HDF5/rain_mass_sfc-<SIM>-AS-2012-01-01-000000-g1.h5"
RmassVname = "/rain"
RdiamFtemp = "HDF5/<SIM>/HDF5/rain_diam_sfc-<SIM>-AS-2012-01-01-000000-g1.h5"
RdiamVname = "/rain_diam"

# Output dataset names
EfluxVname = "/energy_flux"
DomainEfluxVname = "/domain_energy_flux"

DefaultRainDiam = 1.0e-3  # Rain mean diameter = 1 mm when IRAIN = 1 in the RAMSIN file
                          # Set units to meters
Cl = 4186.0  # heat capacity of water, J/kg-water/deg C

VrainCoeff = 144.0  # Power law for rain fall speed, from IPLAWS = 2 tables in RAMS, mic_init.f90
VrainPower = 0.497

DtLong = 5.0        # seconds (main time step)

HiResCellHeight =  50.0   # meters
LoResCellHeight = 101.0   # meters

for isim in range(Nsims):
    Sim = SimList[isim]

    # For RCE_3km_2mom_dm_lrz sim, cell height is twice that of other sims.
    if (re.search("_lrz", Sim)):
        CellHeight = LoResCellHeight
    else:
        CellHeight = HiResCellHeight

    print("************************************************************************************")
    print("Generating rain drop energy for simulation: {0:s}".format(Sim))
    print("  CellHeight: ", CellHeight)
    print('')

    # For single moment, use a default rain drop diameter.
    # For double moment, read the diagnosed rain drop diameter from the REVU output.
    UseDefaultDiam = re.search("1mom", Sim) is not None

    # Form input file names
    TdiffFname = TdiffFtemp.replace("<SIM>", Sim)
    DensFname  = DensFtemp.replace("<SIM>",  Sim)
    RmassFname = RmassFtemp.replace("<SIM>", Sim)
    if not UseDefaultDiam:
        RdiamFname = RdiamFtemp.replace("<SIM>", Sim)

    # Open the input files
    TdiffFile = h5py.File(TdiffFname, mode='r')
    DensFile  = h5py.File(DensFname,  mode='r')
    RmassFile = h5py.File(RmassFname, mode='r')
    if not UseDefaultDiam:
        RdiamFile = h5py.File(RdiamFname, mode='r')

    # Read in the coordinates from the temperature difference file
    print("  Reading: {0:s} ({1:s})".format(TdiffFname, Xname))
    print("  Reading: {0:s} ({1:s})".format(TdiffFname, Yname))
    print("  Reading: {0:s} ({1:s})".format(TdiffFname, Zname))
    print("  Reading: {0:s} ({1:s})".format(TdiffFname, Tname))
    print("")

    X = TdiffFile[Xname][...]
    Y = TdiffFile[Yname][...]
    Z = TdiffFile[Zname][...]
    T = TdiffFile[Tname][...]

    Nx = len(X)
    Ny = len(Y)
    Nz = len(Z)
    Nt = len(T)

    # Inform user that we are reading the input files
    print("  Reading: {0:s} ({1:s})".format(TdiffFname, TdiffVname))
    print("  Reading: {0:s} ({1:s})".format(DensFname,  DensVname))
    print("  Reading: {0:s} ({1:s})".format(RmassFname, RmassVname))
    if not UseDefaultDiam:
        print("  Reading: {0:s} ({1:s})".format(RdiamFname, RdiamVname))
    print("")

    # Open the output file
    OutFname = "DIAGS/rain_drop_energy_{0:s}.h5".format(Sim)
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
    EfluxDset = h5u.DsetCoards(EfluxVname, 3, [Nt, Ny, Nx], chunks=(1, Ny, Nx))
    EfluxDset.Create(Ofile)

    DomainEfluxDset = h5u.DsetCoards(DomainEfluxVname, 1, [Nt])
    DomainEfluxDset.Create(Ofile)

    # Process one time step at a time
    for i in range(Nt):
        Tdiff = np.squeeze(TdiffFile[TdiffVname][i,...])           # deg C
        Dens  = np.squeeze(DensFile[DensVname][i,...])             # kg-air/m^3
        Rmass = np.squeeze(RmassFile[RmassVname][i,...]) * 1.0e-3  # kg-rain/kg-air
        if UseDefaultDiam:
            Rdiam = np.ones((Ny, Nx)) * DefaultRainDiam
        else:
            Rdiam = np.squeeze(RdiamFile[RdiamVname][i,...]) * 1.0e-3 # meters

        # Calculate the portion of rain mass that leaves the grid cell (Alpha). Assume the
        # mass is evenly distributed across the cell. Then, according to rain diameter,
        # calculate how fast the rain is falling. Combine the fall speed of rain with
        # the vertical distance of the grid cell to calculate how much of the rain
        # mass exits the cell at the bottom.
        #
        # Use the RAMS power law formula to calculate the fall speed of rain. This is
        # taken from the tables in mic_init.f90 for IPLAWS = 2.
        Vrain = VrainCoeff * np.power(Rdiam, VrainPower)  # m/s

        # Portion is the distance rain falls divided by the cell height. If the
        # rain falls more than the cell height (Alpha > 1), then reset Alpha to one.
        Alpha = (Vrain * DtLong) / CellHeight
        Alpha[Alpha > 1.0] = 1.0

        # Calculate energy flux in every horizontal grid cell, then take mean across
        # the domain. Include zeros so that the fact that it is not raining everywhere
        # is taken into account.
        #
        # Treat all cells with rain mass < 1e-5 kg/kg (0.01 g/kg) as zero rain mass which
        # means zero contribution to the domain energy flux. Zero out Eflux where Rmass
        # is <= 1.0e-5, and then take domain mean of eflux.
        Eflux = Cl * Rmass * Dens * Alpha * Tdiff * Vrain
        Eflux[(Rmass <= 1.0e-5)] = 0.0
        DomainEflux = np.mean(Eflux)

        Ofile[EfluxVname][i,...] = Eflux
        Ofile[DomainEfluxVname][i] = DomainEflux

        if ((i % 10) == 0):
            print("    Processing time step: {0:d}".format(i))

    print("")

    # close input files
    TdiffFile.close()
    DensFile.close()
    RmassFile.close()
    if not UseDefaultDiam:
        RdiamFile.close()

    # Attach dimensions to output datasets
    print("  Writing {0:s} ({1:s})".format(OutFname, EfluxVname))
    EfluxDset.AttachDims(Ofile, Tdim, Ydim, Xdim)
    print("  Writing {0:s} ({1:s})".format(OutFname, DomainEfluxVname))
    DomainEfluxDset.AttachDims(Ofile, Tdim)

    print("")

    Ofile.close()
