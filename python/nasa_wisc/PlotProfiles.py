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
#    [ [ 'RCE_1km', 'RCE_1km_SM', 'RCE_1km_DM', 'RCE_1km_DP' ], r'$IPLAWS = 0$', 'RceIpl0' ],
#    [ [ 'RCE_1km_IPL2', 'RCE_1km_SM_IPL2', 'RCE_1km_DM_IPL2', 'RCE_1km_DP_IPL2' ], r'$IPLAWS = 2$', 'RceIpl2' ],
    [ [ 'RCE_3km_1mom', 'RCE_3km_2mom', 'RCE_3km_2mom_db', 'RCE_3km_2mom_dm' ], '', 'Rce3km' ],
    ]
Nplots = len(PlotList)

VarList = [
    [ 'DIAGS/profiles_<SIM>.h5', 'avg_theta_prof',       [ 290, 420 ],     r'$\theta\ (K)$',               'Theta'      ],
    [ 'DIAGS/profiles_<SIM>.h5', 'avg_total_cond_prof',  [ -0.001, 0.08 ], r'$Total Cond.\ (g\ kg^{-1})$', 'TotalCond'  ],
    [ 'DIAGS/profiles_<SIM>.h5', 'avg_vapor_prof',       [ -0.5, 20 ],     r'$Vapor\ (g\ kg^{-1})$',       'Vapor'      ],
    ]
Nvars = len(VarList)

Z = np.array(conf.SetZcoords()) / 1000.0 # km
Nz = len(Z)

for iplot in range(Nplots):
    SimList = PlotList[iplot][0]
    Nsims = len(SimList)
    Pname = PlotList[iplot][1]
    Pfprefix = PlotList[iplot][2]

    for ivar in range(Nvars):
        InFileTemplate = VarList[ivar][0]
        Var = VarList[ivar][1]
        Xlims = VarList[ivar][2]
        Xlabel = VarList[ivar][3]
        Oname = VarList[ivar][4]

        print("Generating plot for variable: {0:s}".format(Var))

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
                VARS = np.zeros((Nsims, Nz), dtype=np.float_)

            VARS[isim,...] = InFile[Vname][...]

            InFile.close()

        print("")

        # Make the plot
        Pdir = 'Plots.py'
        Ylims = [ 0, 18 ] # km
        Ylabel = r'$Height\ (km)$'
        
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
        
        plu.PlotLine(Ax, VARS.transpose(), Z, Ptitle, Xaxis, Yaxis, Legend, SimColors)
        
        
        OutFile = "{0:s}/{1:s}Profile{2:s}.png".format(Pdir, Pfprefix, Oname)
        print("Writing: {0:s}".format(OutFile))
        Fig.savefig(OutFile)
        plt.close()

        print("")
        
