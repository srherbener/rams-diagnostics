#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))

import matplotlib.pyplot as plt
import numpy as np
import h5py
from metpy import calc as metc
from metpy import plots as metp
# Access the pint unit registry through metpy so that the metpy
# routines are using the same registry.
from metpy.units import units
import ConfigRce as conf

SimList = [
    "RCE_3km_1mom",
    "RCE_3km_1mom_db",
    "RCE_3km_1mom_db_udef",
    "RCE_3km_1mom_db_rlongup",
    "RCE_3km_1mom_dm",
    "RCE_3km_2mom",
    "RCE_3km_2mom_db",
    "RCE_3km_2mom_db_udef",
    "RCE_3km_2mom_db_rlongup",
    "RCE_3km_2mom_dm",
    "RCE_3km_2mom_dm_lrz",
    ]
Nsims = len(SimList)

InFileTemplate = "DIAGS/profiles_<SIM>.h5"
OutFileTemplate = "Plots.py/RceSkewT_<SIM>.png"

TempVname = "/tempk_prof"
RhVname   = "/relhum_prof"

Zname = "/z_coords"

LabelScheme = conf.SetLabelScheme()

for isim in range(Nsims):
    Sim = SimList[isim]
    print("********** Creating sounding plots for simulation: {0:s}".format(Sim))
    print("")
 
    InFname = InFileTemplate.replace("<SIM>", Sim)
    InFile = h5py.File(InFname, mode='r')

    print("  Reading: {0:s} ({1:s})".format(InFname, TempVname))
    print("  Reading: {0:s} ({1:s})".format(InFname, RhVname))
    print("  Reading: {0:s} ({1:s})".format(InFname, Zname))
    print("")

    # MetPy wants quantities in the pint format which stores an
    # array with the associate units.
    T  = (InFile[TempVname][...] * units('kelvin')).to(units('degC')) # input temp is K, convert to deg C
    RH = (InFile[RhVname][...] / 100.0 * units('radian')).to('')      # convert to dimensionless, use
                                                                      # radian initially since it's a
                                                                      # dimensionless quantity
    Z = InFile[Zname][...] * units('meter')
    Nz = len(Z)

    InFile.close()

    # Convert the heights to pressure values
    P = metc.height_to_pressure_std(Z)  # P comes back in millibars

    # Calculate the dew point temperature
    Td = metc.dewpoint_rh(T, RH)       # Td comes back in deg C

    # Make the plot
    Ptitle = "{0:s} (averaged over final 20 days)".format(LabelScheme[Sim])
    Fig = plt.figure(figsize=(9, 9))

    skewt = metp.SkewT(Fig, rotation=30)

    skewt.plot(P, T, 'r')
    skewt.plot(P, Td, 'g')
    skewt.ax.set_xlim(-80, 30)
    skewt.ax.set_ylim(1000, 50)
    skewt.ax.title.set_text(Ptitle)

    skewt.plot_dry_adiabats()
    skewt.plot_moist_adiabats()
    skewt.plot_mixing_lines()

    OutFname = OutFileTemplate.replace("<SIM>", Sim)
    print("  Writing: {0:s}".format(OutFname))
    Fig.savefig(OutFname)
    plt.close()

    print("")
