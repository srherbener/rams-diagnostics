#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))

import matplotlib.pyplot as plt
import numpy as np
import PlotUtils as plu
import h5py

CaseList = [ [ 'TSD_SAL_DUST',      'SD',   'black'   ],
             [ 'TSD_SAL_NODUST',    'SND',  'cyan'    ],
             [ 'TSD_NONSAL_DUST',   'NSD',  'magenta' ],
             [ 'TSD_NONSAL_NODUST', 'NSND', 'orange'  ] ]
Ncases = len(CaseList)

# Read in data for the Vt max time series panel
LegText = [ ]
Colors = [ ]
for i in range(0, Ncases):
    Case  = CaseList[i][0]
    Label = CaseList[i][1]
    Color = CaseList[i][2]

    InFname = "DIAGS/storm_meas_{0:s}.h5".format(Case)
    InVname = "/max_wind_t"
    print("Reading {0:s} ({1:s})".format(InFname, InVname))

    InFile = h5py.File(InFname, mode='r')

    # grab the time coordinate values, and initialize the wind array
    if (i == 0):
        T = InFile['/t_coords'][:] / 3600 - 42  # sim time in hours starting with zero
        Nt = len(T)
        WindTs = np.zeros( [ Nt, Ncases ], dtype=np.float_)

    WindTs[ :, i ] = InFile[InVname][:]
    InFile.close()

    LegText.append(Label)
    Colors.append(Color)

print("")

# Read in PSAP and SAP factors
InFname = "DIAGS/fs_factors.h5"
InFile = h5py.File(InFname, mode='r')

InVname = "/s_max_wind_t_bar_factors"
print("Reading {0:s} ({1:s})".format(InFname, InVname))
SalWindH1 = InFile[InVname][:,0]
SalWindH2 = InFile[InVname][:,1]

InVname = "/ps_max_wind_t_bar_factors"
print("Reading {0:s} ({1:s})".format(InFname, InVname))
PreSalWindH1 = InFile[InVname][:,0]
PreSalWindH2 = InFile[InVname][:,1]

InFile.close()

print("")

# Create 3 panel plot
#    top half: Vt max time series
#    lower left quarter: PSAP bar graph
#    lower right quarter: SAP bar graph
Pdir = 'Plots.py'

Nbars = SalWindH1.size
Xbars = np.arange(1, Nbars+1) # x values for bar graphs: 1 through Nbars

Fig = plt.figure()

AxVt   = Fig.add_axes([ 0.1, 0.6, 0.7, 0.3 ])
AxPsap = Fig.add_axes([ 0.1, 0.1, 0.3, 0.3 ])
AxSap  = Fig.add_axes([ 0.6, 0.1, 0.3, 0.3 ])

# Vt max time series, top half
Xaxis = plu.AxisConfig('x', [ T[0], T[-1] ], '')
Xaxis.ticks = [ 6, 18, 30, 42, 54 ]
Xaxis.ticklabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n23Aug" ] 
Xaxis.fontsize = 14

Yaxis = plu.AxisConfig('y', [ 5, 23 ], r'$Speed (ms^{-1})$')
Yaxis.fontsize = 14

Ptitle = plu.TitleConfig('a', r'$V_t (max)$')
Ptitle.fontsize = 24

Legend = plu.LegendConfig(LegText, 'upper left')
Legend.bbox = [ 1.02, 1.0 ]
Legend.fontsize = 12
Legend.ncol = 1

plu.PlotLine(AxVt, T, WindTs, Ptitle, Xaxis, Yaxis, Legend, Colors)

# Mark the PRESAL (10-30 h) and SAL (40-60 h) time periods
AxVt.plot([ 10, 30 ], [ 7, 7 ], color='black', linewidth=2);
AxVt.plot([ 10, 10 ], [ 6, 8 ], color='black', linewidth=2);
AxVt.plot([ 30, 30 ], [ 6, 8 ], color='black', linewidth=2);
AxVt.text(15, 8.4, 'PSAP', color='black', fontsize=14);
 
AxVt.plot([ 40, 60 ], [ 7, 7 ], color='black', linewidth=2);
AxVt.plot([ 40, 40 ], [ 6, 8 ], color='black', linewidth=2);
AxVt.plot([ 60, 60 ], [ 6, 8 ], color='black', linewidth=2);
AxVt.text(48, 8.4, 'SAP', color='black', fontsize=14);

# For the bar graphs
BarColors = [ 'orange', 'cyan', 'magenta', 'lightgreen', 'black' ]
BarYlim = [ 13, 18 ]
BarYlabel = r'$Speed (ms^{-1})$'

# PSAP bar graph
Xaxis = plu.AxisConfig('x', [ 0, Nbars+1 ], '')
Xaxis.ticks = range(1, Nbars+1)
Xaxis.ticklabels = [ 'NSND', 'F1', 'F2', 'F12', 'SD' ]
Xaxis.fontsize = 14

Yaxis = plu.AxisConfig('y', BarYlim, BarYlabel)
Yaxis.fontsize = 14

Ptitle = plu.TitleConfig('b', 'PSAP')
Ptitle.fontsize = 24

Legend = plu.LegendConfig([ ], 'none')

plu.PlotSplitBgraph(AxPsap, Xbars, PreSalWindH1, PreSalWindH2, Ptitle, Xaxis, Yaxis, Legend, BarColors)

# SAP bar graph
Xaxis = plu.AxisConfig('x', [ 0, Nbars+1 ], '')
Xaxis.ticks = range(1, Nbars+1)
Xaxis.ticklabels = [ 'NSND', 'F1', 'F2', 'F12', 'SD' ]
Xaxis.fontsize = 14

Yaxis = plu.AxisConfig('y', BarYlim, BarYlabel)
Yaxis.fontsize = 14

Ptitle = plu.TitleConfig('c', 'SAP')
Ptitle.fontsize = 24

Legend = plu.LegendConfig([ ], 'none')

plu.PlotSplitBgraph(AxSap, Xbars, SalWindH1, SalWindH2, Ptitle, Xaxis, Yaxis, Legend, BarColors)




OutFile = "{0:s}/FsFigWindMax.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)

