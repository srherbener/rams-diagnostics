#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import matplotlib.pyplot as plt
import numpy as np
import h5py
import PlotUtils as plu
import ConfigTsd as conf

# color and label schemes for plotting
ColorScheme = conf.SetColorScheme()
LabelScheme = conf.SetLabelScheme()

# Read in factor separation data
# Arrays in files have the time series in rows
# In the sims var:
#       row      sim
#        1      NSND  (S0)
#        2      SND   (S1)
#        3      NSD   (S2)
#        4      SD    (S12)
#
# In the factors var:
#       row     factor
#        1       F0
#        2       F1
#        3       F2
#        4       F12

# The following need to stay in sync with array
# as described above.
SimsList = [ 
    'TSD_NONSAL_NODUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_SAL_DUST'
    ]

FacsList = [ 
    'TSD_NONSAL_NODUST',
    'FAC_SAL',
    'FAC_DUST',
    'FAC_INT',
    'TSD_SAL_DUST'
    ]

InFname = 'DIAGS/fs_factor_time_series.h5'
InFile = h5py.File(InFname, mode='r')

InVname = '/max_wind_t_sims'
print("Reading {0:s} ({1:s})".format(InFname, InVname))
Dset = InFile[InVname]
( Ns, Nt ) = Dset.shape
Sims = np.zeros([ Ns, Nt ])
Sims = Dset[:,:]

InVname = '/max_wind_t_factors'
print("Reading {0:s} ({1:s})".format(InFname, InVname))
Dset = InFile[InVname]
( Nf, Nt ) = Dset.shape
RawFactors = np.zeros([ Nf, Nt ])
RawFactors = Dset[:,:]

InVname = '/t_coords'
print("Reading {0:s} ({1:s})".format(InFname, InVname))
Dset = InFile[InVname]
T = np.zeros([ Nt ])
T = Dset[:] / 3600 - 42 # sim time in hrs

InFile.close()
print("")

# Form an array where each column is a time series ordered:
#   [ NSND F1 F2 F12 SD ]
#
#  NSND is row 1 of Sims
#  F1   is row 2 of Factors
#  F1   is row 3 of Factors
#  F12  is row 4 of Factors
#  SD   is row 4 of Sims

Factors = np.zeros([ Nt, 5 ])
Factors[:,0] = Sims[0,:]
Factors[:,1] = RawFactors[1,:]
Factors[:,2] = RawFactors[2,:]
Factors[:,3] = RawFactors[3,:]
Factors[:,4] = Sims[3,:]


# Create 2 panel plot
#    top half: Vt max time series
#    bottom half: Vt factors
Pdir = 'Plots.py'

SimColors = [ ]
SimLegText = [ ]
for Case in SimsList:
    SimColors.append(ColorScheme[Case])
    SimLegText.append(LabelScheme[Case])

FacColors = [ ]
FacLegText = [ ]
for Case in FacsList:
    FacColors.append(ColorScheme[Case])
    FacLegText.append(LabelScheme[Case])


Fig = plt.figure()

AxSims    = Fig.add_axes([ 0.1, 0.6, 0.7, 0.3 ])
AxFactors = Fig.add_axes([ 0.1, 0.1, 0.7, 0.3 ])

# Vt max time series, top half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ] 
Xaxis.fontsize = 14

Yaxis = plu.AxisConfig('y', [ 5, 23 ], r'$Speed (ms^{-1})$')
Yaxis.fontsize = 14

Ptitle = plu.TitleConfig('a', r'$V_t (max): Sims$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(SimLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

# Sims needs to be transposed for the plotting routine
plu.PlotLine(AxSims, T, Sims.transpose(), Ptitle, Xaxis, Yaxis, Legend, SimColors)

# Mark the PRESAL (10-30 h) and SAL (40-60 h) time periods
AxSims.plot([ 10, 30 ], [ 7, 7 ], color='black', linewidth=2);
AxSims.plot([ 10, 10 ], [ 6, 8 ], color='black', linewidth=2);
AxSims.plot([ 30, 30 ], [ 6, 8 ], color='black', linewidth=2);
AxSims.text(15, 8.4, 'PSAP', color='black', fontsize=14);
 
AxSims.plot([ 40, 60 ], [ 7, 7 ], color='black', linewidth=2);
AxSims.plot([ 40, 40 ], [ 6, 8 ], color='black', linewidth=2);
AxSims.plot([ 60, 60 ], [ 6, 8 ], color='black', linewidth=2);
AxSims.text(48, 8.4, 'SAP', color='black', fontsize=14);



# Factor time series, bottom half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ] 
Xaxis.fontsize = 14

Yaxis = plu.AxisConfig('y', [ -5, 23 ], r'$Speed (ms^{-1})$')
Yaxis.fontsize = 14

Ptitle = plu.TitleConfig('b', r'$V_t (max): Factors$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(FacLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxFactors, T, Factors, Ptitle, Xaxis, Yaxis, Legend, FacColors)

# Add y = 0 reference
AxFactors.axhline(y=0, color='black', linestyle='--')


OutFile = "{0:s}/FsFigWindMaxTseries.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)

