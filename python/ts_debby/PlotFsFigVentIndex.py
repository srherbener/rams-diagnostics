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

InFnameTmpl = 'DIAGS/vent_index_<SIM>.h5'
WshrVname = '/sm_wind_shear'
EdefVname = '/sm_entropy_deficit'
PintVname = '/sm_pot_intensity'
VindVname = '/sm_vent_index'
Tvname = '/t_coords'

SimLegText = []
SimColors = []
for i in range(Nsims):
    Sim = SimsList[i]
    Label = LabelScheme[Sim]
    Color = ColorScheme[Sim]

    InFname = InFnameTmpl.replace("<SIM>", Sim)
    InFile = h5py.File(InFname, mode='r')

    print("Reading: {0:s} ({1:s})".format(InFname, VindVname))
    print("Reading: {0:s} ({1:s})".format(InFname, WshrVname))
    print("Reading: {0:s} ({1:s})".format(InFname, EdefVname))
    print("Reading: {0:s} ({1:s})".format(InFname, PintVname))
    VindVar = InFile[VindVname][...]
    WshrVar = InFile[WshrVname][...]
    EdefVar = InFile[EdefVname][...]
    PintVar = InFile[PintVname][...]

    if (i == 0):
        print("  Reading: {0:s} ({1:s})".format(InFname, Tvname))
        T = InFile[Tvname][...] / 3600 - 42 # Convert seconds to sim time
        Nt = len(T)

        # Allocate an array to hold multiple lines
        WshrSims = np.zeros((Nt, Nsims), dtype=np.float_)
        EdefSims = np.zeros((Nt, Nsims), dtype=np.float_)
        PintSims = np.zeros((Nt, Nsims), dtype=np.float_)
        VindSims = np.zeros((Nt, Nsims), dtype=np.float_)

    WshrSims[:,i] = WshrVar
    EdefSims[:,i] = EdefVar
    PintSims[:,i] = PintVar
    VindSims[:,i] = VindVar
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
WshrFactors = np.zeros((Nt, Nfacs), dtype=np.float_)
EdefFactors = np.zeros((Nt, Nfacs), dtype=np.float_)
PintFactors = np.zeros((Nt, Nfacs), dtype=np.float_)
VindFactors = np.zeros((Nt, Nfacs), dtype=np.float_)
for i in range(Nfacs):
    Factor = FacsList[i]
    Label = LabelScheme[Factor]
    Color = ColorScheme[Factor]

    if (Factor == 'TSD_NONSAL_NODUST'):
      WshrFactors[:,i] = WshrSims[:,0]
      EdefFactors[:,i] = EdefSims[:,0]
      PintFactors[:,i] = PintSims[:,0]
      VindFactors[:,i] = VindSims[:,0]
    elif (Factor == 'FAC_SAL'):
      WshrFactors[:,i] = WshrSims[:,1] - WshrSims[:,0]
      EdefFactors[:,i] = EdefSims[:,1] - EdefSims[:,0]
      PintFactors[:,i] = PintSims[:,1] - PintSims[:,0]
      VindFactors[:,i] = VindSims[:,1] - VindSims[:,0]
    elif (Factor == 'FAC_DUST'):
      WshrFactors[:,i] = WshrSims[:,2] - WshrSims[:,0]
      EdefFactors[:,i] = EdefSims[:,2] - EdefSims[:,0]
      PintFactors[:,i] = PintSims[:,2] - PintSims[:,0]
      VindFactors[:,i] = VindSims[:,2] - VindSims[:,0]
    elif (Factor == 'FAC_INT'):
      WshrFactors[:,i] = WshrSims[:,3] - (WshrSims[:,1] + WshrSims[:,2]) + WshrSims[:,0]
      EdefFactors[:,i] = EdefSims[:,3] - (EdefSims[:,1] + EdefSims[:,2]) + EdefSims[:,0]
      PintFactors[:,i] = PintSims[:,3] - (PintSims[:,1] + PintSims[:,2]) + PintSims[:,0]
      VindFactors[:,i] = VindSims[:,3] - (VindSims[:,1] + VindSims[:,2]) + VindSims[:,0]
    elif (Factor == 'TSD_SAL_DUST'):
      WshrFactors[:,i] = WshrSims[:,3]
      EdefFactors[:,i] = EdefSims[:,3]
      PintFactors[:,i] = PintSims[:,3]
      VindFactors[:,i] = VindSims[:,3]

    FacLegText.append(Label)
    FacColors.append(Color)


# Wind shear plot
#    top half: Sims
#    bottom half: Factors
Pdir = 'Plots.py'
YlimSims = [ 0, 8 ]
YlimFacs = [ -2, 8 ]

Fig = plt.figure()

AxSims    = Fig.add_axes([ 0.1, 0.6, 0.7, 0.3 ])
AxFactors = Fig.add_axes([ 0.1, 0.1, 0.7, 0.3 ])

# Simulations, top half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', YlimSims, r'Shear')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('a', r'$PSAP:\ Wind\ Shear:\ Sims$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(SimLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxSims, T, WshrSims, Ptitle, Xaxis, Yaxis, Legend, SimColors)


# Factors, bottom half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', YlimFacs, r'Shear')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('b', r'$PSAP:\ Wind\ Shear:\ Factors$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(FacLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxFactors, T, WshrFactors, Ptitle, Xaxis, Yaxis, Legend, FacColors)

# Add y = 0 reference
AxFactors.axhline(y=0, color='black', linestyle='--')


OutFile = "{0:s}/FsFigWindShear.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)
plt.close()

# Entropy deficit plot
#    top half: Sims
#    bottom half: Factors
Pdir = 'Plots.py'
YlimSims = [ 0, 1.8 ]
YlimFacs = [ -0.2, 1.8 ]

Fig = plt.figure()

AxSims    = Fig.add_axes([ 0.1, 0.6, 0.7, 0.3 ])
AxFactors = Fig.add_axes([ 0.1, 0.1, 0.7, 0.3 ])

# Simulations, top half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', YlimSims, r'ED')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('a', r'$PSAP:\ Entropy\ Def.:\ Sims$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(SimLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxSims, T, EdefSims, Ptitle, Xaxis, Yaxis, Legend, SimColors)


# Factors, bottom half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', YlimFacs, r'ED')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('b', r'$PSAP:\ Entropy\ Def.:\ Factors$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(FacLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxFactors, T, EdefFactors, Ptitle, Xaxis, Yaxis, Legend, FacColors)

# Add y = 0 reference
AxFactors.axhline(y=0, color='black', linestyle='--')


OutFile = "{0:s}/FsFigEntropyDeficit.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)
plt.close()

# Potential intensity plot
#    top half: Sims
#    bottom half: Factors
Pdir = 'Plots.py'
YlimSims = [ 40, 60 ]
YlimFacs = [ -2, 2 ]

Fig = plt.figure()

AxSims    = Fig.add_axes([ 0.1, 0.6, 0.7, 0.3 ])
AxFactors = Fig.add_axes([ 0.1, 0.1, 0.7, 0.3 ])

# Simulations, top half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', YlimSims, r'PI')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('a', r'$PSAP:\ Pot.\ Intensity:\ Sims$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(SimLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxSims, T, PintSims, Ptitle, Xaxis, Yaxis, Legend, SimColors)


# Factors, bottom half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', YlimFacs, r'PI')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('b', r'$PSAP:\ Pot.\ Intensity:\ Factors$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(FacLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxFactors, T, PintFactors, Ptitle, Xaxis, Yaxis, Legend, FacColors)

# Add y = 0 reference
AxFactors.axhline(y=0, color='black', linestyle='--')


OutFile = "{0:s}/FsFigPotentialIntensity.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)
plt.close()


# Ventilation Index plot
#    top half: Sims
#    bottom half: Factors
Pdir = 'Plots.py'
YlimSims = [ 0, 0.21 ]
YlimFacs = [ -0.05, 0.21 ]

Fig = plt.figure()

AxSims    = Fig.add_axes([ 0.1, 0.6, 0.7, 0.3 ])
AxFactors = Fig.add_axes([ 0.1, 0.1, 0.7, 0.3 ])

# Simulations, top half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', YlimSims, r'VI')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('a', r'$PSAP:\ Vent.\ Index:\ Sims$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(SimLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxSims, T, VindSims, Ptitle, Xaxis, Yaxis, Legend, SimColors)


# Factors, bottom half
Xaxis = plu.AxisConfig('x', [ 0, 62 ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ]
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', YlimFacs, r'VI')
Yaxis.fontsize = 12

Ptitle = plu.TitleConfig('b', r'$PSAP:\ Vent.\ Index:\ Factors$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(FacLegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxFactors, T, VindFactors, Ptitle, Xaxis, Yaxis, Legend, FacColors)

# Add y = 0 reference
AxFactors.axhline(y=0, color='black', linestyle='--')


OutFile = "{0:s}/FsFigVentilationIndex.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)
plt.close()

