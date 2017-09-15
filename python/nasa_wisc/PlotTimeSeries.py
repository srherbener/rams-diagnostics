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

    [ [ 'RCE_3km_1mom', 'RCE_3km_1mom_db', 'RCE_3km_1mom_dm' ], '', 'Rce3km1mom', 0, 120 ],
    [ [ 'RCE_3km_2mom_db', 'RCE_3km_2mom_dm', 'RCE_3km_2mom_dm_lrz' ], '', 'Rce3km2mom', 0, 120 ],
    ]
Nplots = len(PlotList)

VarList = [
    [ 'DIAGS/dom_avg_pcprr_<SIM>.h5',      'pcprr',      [   0,  12 ], r'$PR\ (mm\ day^{-1})$', 'Prate'  ],
    [ 'DIAGS/dom_avg_vint_vapor_<SIM>.h5', 'vint_vapor', [   0,  60 ], r'$PW\ (mm)$',           'Pwater' ],

    [ 'DIAGS/dom_avg_sfc_lat_<SIM>.h5',  'sfc_lat',    [   0, 140 ], r'$LHF\ (W\ m^{2})$',    'Lhf'    ],
    [ 'DIAGS/dom_avg_sfc_sens_<SIM>.h5', 'sfc_sens',   [   0,  20 ], r'$SHF\ (W\ m^{2})$',    'Shf'    ],

    [ 'DIAGS/dom_avg_top_lwup_<SIM>.h5', 'top_lwup',   [ 150, 350 ], r'$OLR\ (W\ m^{2})$',    'Olr'    ],
    [ 'DIAGS/dom_avg_top_swdn_<SIM>.h5', 'top_swdn',   [ 400, 420 ], r'$INSOLATION\ (W\ m^{2})$',    'Insol'    ],

    [ 'DIAGS/eq_meas_<SIM>.h5', 'rad_flux_div',      [   0,  150 ], r'$QRAD\ (W\ m^{2})$',   'Qrad'   ],
    [ 'DIAGS/eq_meas_<SIM>.h5', 'therm_heat_flux',   [   0,  150 ], r'$THF\ (W\ m^{2})$',    'Thf'    ],
    ]
Nvars = len(VarList)

Tvname = '/t_coords'

# plot config
Pdir = 'Plots.py'
Xlabel = r'$Simulation\ Time\ (days)$'

for iplot in range(Nplots):
    SimList  = PlotList[iplot][0]
    Nsims = len(SimList)
    Pname    = PlotList[iplot][1]
    Pfprefix = PlotList[iplot][2]
    Tstart   = PlotList[iplot][3]
    Tend     = PlotList[iplot][4]

    # plot config
    Xlims = [ Tstart, Tend ]
    T1 = Tstart * 24  # time steps are in hours
    T2 = Tend * 24

    for ivar in range(Nvars):
        InFileTemplate = VarList[ivar][0]
        Var    = VarList[ivar][1]
        Ylims  = VarList[ivar][2]
        Ylabel = VarList[ivar][3]
        Oname  = VarList[ivar][4]

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
                T = np.arange(T1, T2+1) / 24 # convert hours to days
                Nt = len(T)
                VARS = np.ones((Nsims, Nt), dtype=np.float_) * np.nan

            IN_VAR = InFile[Vname][T1:T2+1,...]
            VarNt = len(IN_VAR)
            if (Vname == "/pcprr"):
                IN_VAR = IN_VAR * 24 # convert mm/h to mm/day

            VARS[isim,0:VarNt] = IN_VAR

            InFile.close()

        print("")

        # Make the plot
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
        Legend.fontsize = 9
        Legend.ncol = 1
        
        plu.PlotLine(Ax, T, VARS.transpose(), Ptitle, Xaxis, Yaxis, Legend, SimColors)
        
        
        OutFile = "{0:s}/{1:s}TimeSeries{2:s}.png".format(Pdir, Pfprefix, Oname)
        print("Writing: {0:s}".format(OutFile))
        Fig.savefig(OutFile)
        plt.close()

        print("")
        
