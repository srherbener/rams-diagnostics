#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import PlotTcXsection as ptx


InFname = 'DIAGS/storm_xsections_<CASE>.h5'

## Max Vt
#SimCspecs = [ 0, 20, 21 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
#FactCspecs = [ -4, 4, 21 ]
#VtPsap = ptx.StormXsection(InFname, '/all_ps_speed_t', 'V_t', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigVtPsapFactors.png')
#VtPsap.PlotXsection()
#VtSap = ptx.StormXsection(InFname, '/all_s_speed_t', 'V_t', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigVtSapFactors.png')
#VtSap.PlotXsection()
#
## Ice deposition
#SimCspecs = [ 0, 1.5, 16 ]
#FactCspecs = [ -1, 1, 11 ]
#Ylim = [ 0, 15 ]
#Yticks = [ 0, 3, 6, 9, 12, 15 ]
#IceDepPsap = ptx.StormXsection(InFname, '/all_ps_ice_dep', 'Ice Dep', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigIceDepPsapFactors.png')
#IceDepPsap.Ylim = Ylim
#IceDepPsap.Yticks = Yticks
#IceDepPsap.PlotXsection()
#IceDepSap = ptx.StormXsection(InFname, '/all_s_ice_dep', 'Ice Dep', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigIceDepSapFactors.png')
#IceDepSap.Ylim = Ylim
#IceDepSap.Yticks = Yticks
#IceDepSap.PlotXsection()
#
## rime
#SimCspecs = [ 0, 1, 11 ]
#FactCspecs = [ -1, 1, 11 ]
#Ylim = [ 0, 12 ]
#Yticks =  [ 0, 3, 6, 9, 12 ]
#RimePsap = ptx.StormXsection(InFname, '/all_ps_cloud_rime', 'Rime', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigRimePsapFactors.png')
#RimePsap.Ylim = Ylim
#RimePsap.Yticks = Yticks
#RimePsap.PlotXsection()
#RimeSap = ptx.StormXsection(InFname, '/all_s_cloud_rime', 'Rime', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigRimeSapFactors.png')
#RimeSap.Ylim = Ylim
#RimeSap.Yticks = Yticks
#RimeSap.PlotXsection()

# rain to ice
SimCspecs = [ 0, 0.3, 11 ]
FactCspecs = [ -0.15, 0.15, 11 ]
Ylim = [ 0, 12 ]
Yticks =  [ 0, 3, 6, 9, 12 ]
R2icePsap = ptx.StormXsection(InFname, '/all_ps_rain2ice', 'Rain2Ice', 'PSAP', 
    SimCspecs, FactCspecs, 'Plots.py/FsFigR2icePsapFactors.png')
R2icePsap.Ylim = Ylim
R2icePsap.Yticks = Yticks
R2icePsap.PlotXsection()
R2iceSap = ptx.StormXsection(InFname, '/all_s_rain2ice', 'Rain2Ice', 'SAP', 
    SimCspecs, FactCspecs, 'Plots.py/FsFigR2iceSapFactors.png')
R2iceSap.Ylim = Ylim
R2iceSap.Yticks = Yticks
R2iceSap.PlotXsection()

