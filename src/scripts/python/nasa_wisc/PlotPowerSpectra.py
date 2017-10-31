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
#    [ [ 'RCE_1km', 'RCE_1km_SM', 'RCE_1km_DM', 'RCE_1km_DP' ], r'$IPLAWS = 0$', 'RceIpl0', 0, 50 ],
#    [ [ 'RCE_1km_IPL2', 'RCE_1km_SM_IPL2', 'RCE_1km_DM_IPL2', 'RCE_1km_DP_IPL2' ], r'$IPLAWS = 2$', 'RceIpl2', 0, 50 ],
#    [ [ 'RCE_3km_1mom', 'RCE_3km_2mom', 'RCE_3km_2mom_db', 'RCE_3km_2mom_dm' ], '', 'Rce3km', 0, 80 ],

    [ [ 'RCE_3km_1mom', 'RCE_3km_1mom_db', 'RCE_3km_1mom_db_udef', 'RCE_3km_1mom_db_rlongup', 'RCE_3km_1mom_dm' ], '', 'Rce3km1mom', [ 1.5, 300 ] ],
    [ [ 'RCE_3km_2mom_db', 'RCE_3km_2mom_db_udef', 'RCE_3km_2mom_db_rlongup', 'RCE_3km_2mom_dm', 'RCE_3km_2mom_dm_lrz' ], '', 'Rce3km2mom', [ 1.5, 300 ] ],
    ]
Nplots = len(PlotList)

VarList = [
    [ 'DIAGS/power_spectra_<SIM>.h5',   '/pcprr',      [   1e-3,  1e3 ], r'$PSD(PR^2)$',  'Prate'  ],
    [ 'DIAGS/power_spectra_<SIM>.h5',   '/vint_vapor', [   1e-5,  1e3 ], r'$PSD(PW^2)$',  'Pwater' ],
    [ 'DIAGS/power_spectra_<SIM>.h5',   '/top_lwup',   [   1e-3,  1e4 ], r'$PSD(OLR^2)$', 'Olr'    ],
    ]
Nvars = len(VarList)

# plot config
Pdir = 'Plots.py'
Xlabel = r'$Period (h)$'

for iplot in range(Nplots):
    SimList  = PlotList[iplot][0]
    Pname    = PlotList[iplot][1]
    Pfprefix = PlotList[iplot][2]
    Xlims    = PlotList[iplot][3]

    Nsims = len(SimList)

    for ivar in range(Nvars):
        InFileTemplate = VarList[ivar][0]
        Vname          = VarList[ivar][1]
        Ylims          = VarList[ivar][2]
        Ylabel         = VarList[ivar][3]
        Oname          = VarList[ivar][4]

        print("Generating plot for varialbe: {0:s}".format(Vname))
        PsdVname  = "{0:s}_psd".format(Vname)
        FreqVname = "{0:s}_freq".format(Vname)

        SimLegText = []
        SimColors = []
        for isim in range(Nsims):
            Sim = SimList[isim]
            SimLegText.append(LabelScheme[Sim])
            SimColors.append(ColorScheme[Sim])
         
            InFname = InFileTemplate.replace("<SIM>", Sim)
            InFile = h5py.File(InFname, mode='r')

            print("  Reading: {0:s} ({1:s})".format(InFname, PsdVname))
            PSD = InFile[PsdVname][...]
            if (isim == 0):
                print("  Reading: {0:s} ({1:s})".format(InFname, FreqVname))
                FREQ = InFile[FreqVname][...]
                PER = 1 / FREQ
                Nx = len(PER)

                VARS = np.zeros((Nsims,Nx))

            VARS[isim,...] = PSD

            InFile.close()

        print("")

        # Make the plot
        Fig = plt.figure()
        Ax  = Fig.add_axes([ 0.1, 0.1, 0.7, 0.8 ])
        
        Xaxis = plu.AxisConfig('x', Xlims, Xlabel)
        Xaxis.fontsize = 12
        Xaxis.scale = 'log'
        
        Yaxis = plu.AxisConfig('y', Ylims, Ylabel)
        Yaxis.fontsize = 12
        Yaxis.scale = 'log'
        
        Ptitle = plu.TitleConfig('', Pname)
        Ptitle.fontsize = 24
        
        Legend = plu.LegendConfig(SimLegText, 'upper left')
        Legend.bbox = [ 1.02, 1.0 ]
        Legend.fontsize = 9
        Legend.ncol = 1
        
        plu.PlotLine(Ax, PER, VARS.transpose(), Ptitle, Xaxis, Yaxis, Legend, SimColors)
        
        
        OutFile = "{0:s}/{1:s}PowerSpectrum{2:s}.png".format(Pdir, Pfprefix, Oname)
        print("Writing: {0:s}".format(OutFile))
        Fig.savefig(OutFile)
        plt.close()

        print("")
        
