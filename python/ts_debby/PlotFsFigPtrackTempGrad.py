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
Nsims = len(SimsList)

FacsList = [ 
    'TSD_NONSAL_NODUST',
    'FAC_SAL',
    'FAC_DUST',
    'FAC_INT',
    'TSD_SAL_DUST'
    ]
Nfacs = len(FacsList)

InFnameTmpl = 'DIAGS/ptrack_tgrads_<SIM>.h5'
InVname = '/ps_tempc_nv_p_bar_hgrad_smooth'
Xvname = '/x_coords'

SimLegText = []
SimColors = []
for i in range(Nsims):
    Sim = SimsList[i]
    Label = LabelScheme[Sim]
    Color = ColorScheme[Sim]

    InFname = InFnameTmpl.replace("<SIM>", Sim)
    InFile = h5py.File(InFname, mode='r')

    print("Reading: {0:s} ({1:s})".format(InFname, InVname))
    InVar = InFile[InVname][...] * 1.0e3 # convert to 1e-6 K/m

    if (i == 0):
        print("  Reading: {0:s} ({1:s})".format(InFname, Xvname))
        X = InFile[Xvname][...]
        Nx = len(X)

        # Allocate an array to hold multiple lines
        Sims = np.zeros((Nx, Nsims), dtype=np.float_)

    Sims[:,i] = InVar
    InFile.close()

    SimLegText.append(Label)
    SimColors.append(Color)

print("")

# Form factors for the temp gradient lines. The indices in Sims
# are controlled by the order that SimsList is defined above and needs
# to remain:
#
#   Index      Simulation
#     0          NSND
#     1           SND
#     2           NSD
#     3            SD
#
# The factor order is controlled by FacsList and are formed by:
#
#    Index     Factor            Formula
#      0          NSND      NSND
#      1       FAC_SAL      SND - NSND
#      2      FAC_DUST      NSD - NSND
#      3       FAC_INT      SD - (SND+NSD) + NSND
#      4            SD      SD
#
FacLegText = []
FacColors = []
Factors = np.zeros((Nx, Nfacs), dtype=np.float_)
for i in range(Nfacs):
    Factor = FacsList[i]
    Label = LabelScheme[Factor]
    Color = ColorScheme[Factor]

    if (Factor == 'TSD_NONSAL_NODUST'):
      Factors[:,i] = Sims[:,0]
    elif (Factor == 'FAC_SAL'):
      Factors[:,i] = Sims[:,1] - Sims[:,0]
    elif (Factor == 'FAC_DUST'):
      Factors[:,i] = Sims[:,2] - Sims[:,0]
    elif (Factor == 'FAC_INT'):
      Factors[:,i] = Sims[:,3] - (Sims[:,1] + Sims[:,2]) + Sims[:,0]
    elif (Factor == 'TSD_SAL_DUST'):
      Factors[:,i] = Sims[:,3]

    FacLegText.append(Label)
    FacColors.append(Color)


# Create 2 panel plot
#    top half: Tgrad
#    bottom half: Tgrad factors
Pdir = 'Plots.py'

Fig = plt.figure()

AxSims    = Fig.add_axes([ 0.1, 0.6, 0.7, 0.3 ])
AxFactors = Fig.add_axes([ 0.1, 0.1, 0.7, 0.3 ])

# Vt max time series, top half
Xaxis = plu.AxisConfig('x', [ 0, 1800 ], 'Linear Distance (km)')
Xaxis.ticks = [ 500, 1000, 1500 ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', [ -2.5, 1.5 ], r'${\nabla}\overline{T}\ (10^{-6}\ K\ m^{-1})$')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('a', r'$PSAP:\ {\nabla}\overline{T}:\ Sims$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(SimLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

# Sims needs to be transposed for the plotting routine
plu.PlotLine(AxSims, X, Sims, Ptitle, Xaxis, Yaxis, Legend, SimColors)

# Mark the endpoints of PTRACK
AxSims.text(0, -3, 'C', color='blue', fontsize=16);
AxSims.text(1750, -3, 'D', color='blue', fontsize=16);



# Factor time series, bottom half
Xaxis = plu.AxisConfig('x', [ 0, 1800 ], 'Linear Distance (km)')
Xaxis.ticks = [ 500, 1000, 1500 ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', [ -2.5, 1.5 ], r'${\nabla}\overline{T}\ (10^{-6}\ K\ m^{-1})$')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('b', r'$PSAP:\ {\nabla}\overline{T}:\ Factors$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(FacLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxFactors, X, Factors, Ptitle, Xaxis, Yaxis, Legend, FacColors)

# Mark the endpoints of PTRACK
AxFactors.text(0, -3, 'C', color='blue', fontsize=16);
AxFactors.text(1750, -3, 'D', color='blue', fontsize=16);

# Add y = 0 reference
AxFactors.axhline(y=0, color='black', linestyle='--')


OutFile = "{0:s}/FsFigPtrackTempGrad.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)

