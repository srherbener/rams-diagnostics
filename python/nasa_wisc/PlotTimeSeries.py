#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/nasa_wisc".format(os.environ['HOME']))

import matplotlib.pyplot as plt
import numpy as np
import h5py
import PlotUtils as plu
import ConfigRce as conf

# color and label schemes for plotting
ColorScheme = conf.SetColorScheme()
LabelScheme = conf.SetLabelScheme()

PlotList = [ 
    [ [ 'RCE_1km', 'RCE_1km_SM', 'RCE_1km_DM', 'RCE_1km_DP' ], r'$IPLAWS = 0$', 'RceIpl0' ],
    [ [ 'RCE_1km_IPL2', 'RCE_1km_SM_IPL2', 'RCE_1km_DM_IPL2', 'RCE_1km_DP_IPL2' ], r'$IPLAWS = 2$', 'RceIpl2' ],
    ]
Nplots = len(PlotList)

VarList = [
    [ 'DIAGS/dom_averages_<SIM>.h5', 'avg_olr',       [ 150, 350 ], r'$OLR\ (W\ m^{2})$',    'Olr'    ],
    [ 'DIAGS/dom_averages_<SIM>.h5', 'avg_pcp_rate',  [   0,   8 ], r'$PR\ (mm\ day^{-1})$', 'Prate'  ],
    [ 'DIAGS/dom_averages_<SIM>.h5', 'avg_pcp_water', [  20,  60 ], r'$PW\ (mm)$',           'Pwater' ],
    [ 'DIAGS/dom_averages_<SIM>.h5', 'avg_sfc_lhf',   [   0, 100 ], r'$LHF\ (W\ m^{2})$',    'Lhf'    ],
    [ 'DIAGS/dom_averages_<SIM>.h5', 'avg_sfc_shf',   [   0,  15 ], r'$SHF\ (W\ m^{2})$',    'Shf'    ],
    ]
Nvars = len(VarList)

Tvname = '/t_coords'

for iplot in range(Nplots):
    SimList = PlotList[iplot][0]
    Nsims = len(SimList)
    Pname = PlotList[iplot][1]
    Pfprefix = PlotList[iplot][2]

    for ivar in range(Nvars):
        InFileTemplate = VarList[ivar][0]
        Var = VarList[ivar][1]
        Ylims = VarList[ivar][2]
        Ylabel = VarList[ivar][3]
        Oname = VarList[ivar][4]

        print("Generating plot for varialbe: {0:s}".format(Var))

        Vname = "/{0:s}".format(Var)

        SimLegText = []
        SimColors = []
        for isim in range(Nsims):
            Sim = SimList[isim]
            SimLegText.append(LabelScheme[Sim])
            SimColors.append(ColorScheme[Sim])
         
            InFname = InFileTemplate.replace("<SIM>", Sim)
            print("  Reading: {0:s} ({1:s})".format(InFname, Vname))
            InFile = h5py.File(InFname, mode='r')

            if (isim == 0):
                print("    Reading: {0:s} ({1:s})".format(InFname, Tvname))
                T = InFile[Tvname][...] / 86400.0 # convert sec to day
                Nt = len(T)

                VARS = np.zeros((Nsims, Nt), dtype=np.float_)

            VARS[isim,...] = InFile[Vname][...]

            InFile.close()

        print("")

        # Make the plot
        Pdir = 'Plots.py'
        Xlims = [ 0, 50 ] # 50 days
        Xlabel = r'$Simulation\ Time\ (days)$'
        
        Fig = plt.figure()
        Ax  = Fig.add_axes([ 0.1, 0.1, 0.7, 0.8 ])
        
        Xaxis = plu.AxisConfig('x', Xlims, Xlabel)
        Xaxis.fontsize = 12
        
        Yaxis = plu.AxisConfig('y', Ylims, Ylabel)
        Yaxis.fontsize = 12
        
        Ptitle = plu.TitleConfig('', Pname)
        Ptitle.fontsize = 24
        
        Legend = plu.LegendConfig(SimLegText, 'upper left')
        Legend.bbox = [ 1.02, 1.0 ]
        Legend.fontsize = 12
        Legend.ncol = 1
        
        plu.PlotLine(Ax, T, np.transpose(VARS), Ptitle, Xaxis, Yaxis, Legend, SimColors)
        
        
        OutFile = "{0:s}/{1:s}TimeSeries{2:s}.png".format(Pdir, Pfprefix, Oname)
        print("Writing: {0:s}".format(OutFile))
        Fig.savefig(OutFile)
        plt.close()

        print("")
        
