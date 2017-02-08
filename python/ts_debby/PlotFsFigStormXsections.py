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
#
## rain to ice
#SimCspecs = [ 0, 0.3, 11 ]
#FactCspecs = [ -0.15, 0.15, 11 ]
#Ylim = [ 0, 12 ]
#Yticks =  [ 0, 3, 6, 9, 12 ]
#R2icePsap = ptx.StormXsection(InFname, '/all_ps_rain2ice', 'Rain2Ice', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigR2icePsapFactors.png')
#R2icePsap.Ylim = Ylim
#R2icePsap.Yticks = Yticks
#R2icePsap.PlotXsection()
#R2iceSap = ptx.StormXsection(InFname, '/all_s_rain2ice', 'Rain2Ice', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigR2iceSapFactors.png')
#R2iceSap.Ylim = Ylim
#R2iceSap.Yticks = Yticks
#R2iceSap.PlotXsection()
#
## Theta-E
#SimCspecs = [ 340, 355, 11 ]
#FactCspecs = [ -6, 6, 11 ]
#Ylim = [ 0, 5 ]
#Yticks =  [ 0, 1, 2, 3, 4, 5 ]
#ThetaePsap = ptx.StormXsection(InFname, '/all_ps_theta_e', 'Theta-E', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigThetaePsapFactors.png')
#ThetaePsap.Ylim = Ylim
#ThetaePsap.Yticks = Yticks
#ThetaePsap.PlotXsection()
#ThetaeSap = ptx.StormXsection(InFname, '/all_s_theta_e', 'Theta-E', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigThetaeSapFactors.png')
#ThetaeSap.Ylim = Ylim
#ThetaeSap.Yticks = Yticks
#ThetaeSap.PlotXsection()

# Theta
SimCspecs = [ 290, 340, 11 ]
FactCspecs = [ -3, 3, 11 ]
Ylim = [ 0, 8 ]
Yticks =  [ 0, 2, 4, 6, 8 ]
ThetaPsap = ptx.StormXsection(InFname, '/all_ps_theta', 'Theta', 'PSAP', 
    SimCspecs, FactCspecs, 'Plots.py/FsFigThetaPsapFactors.png')
ThetaPsap.Ylim = Ylim
ThetaPsap.Yticks = Yticks
ThetaPsap.PlotXsection()
ThetaSap = ptx.StormXsection(InFname, '/all_s_theta', 'Theta', 'SAP', 
    SimCspecs, FactCspecs, 'Plots.py/FsFigThetaSapFactors.png')
ThetaSap.Ylim = Ylim
ThetaSap.Yticks = Yticks
ThetaSap.PlotXsection()

## Cloud Evap
#SimCspecs = [ -2, 0, 11 ]
#FactCspecs = [ -0.5, 0.5, 11 ]
#Ylim = [ 0, 5 ]
#Yticks =  [ 0, 1, 2, 3, 4, 5 ]
#CloudEvapPsap = ptx.StormXsection(InFname, '/all_ps_cloud_evap', 'Cloud Evap', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigCloudEvapPsapFactors.png')
#CloudEvapPsap.Ylim = Ylim
#CloudEvapPsap.Yticks = Yticks
#CloudEvapPsap.PlotXsection()
#CloudEvapSap = ptx.StormXsection(InFname, '/all_s_cloud_evap', 'Cloud Evap', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigCloudEvapSapFactors.png')
#CloudEvapSap.Ylim = Ylim
#CloudEvapSap.Yticks = Yticks
#CloudEvapSap.PlotXsection()
#
## Rain Evap
#SimCspecs = [ -1, 0, 11 ]
#FactCspecs = [ -0.5, 0.5, 11 ]
#Ylim = [ 0, 5 ]
#Yticks =  [ 0, 1, 2, 3, 4, 5 ]
#RainEvapPsap = ptx.StormXsection(InFname, '/all_ps_rain_evap', 'Rain Evap', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigRainEvapPsapFactors.png')
#RainEvapPsap.Ylim = Ylim
#RainEvapPsap.Yticks = Yticks
#RainEvapPsap.PlotXsection()
#RainEvapSap = ptx.StormXsection(InFname, '/all_s_rain_evap', 'Rain Evap', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigRainEvapSapFactors.png')
#RainEvapSap.Ylim = Ylim
#RainEvapSap.Yticks = Yticks
#RainEvapSap.PlotXsection()

# Vapor
SimCspecs = [ 0, 20, 11 ]
FactCspecs = [ -2, 2, 11 ]
Ylim = [ 0, 8 ]
Yticks =  [ 0, 2, 4, 6, 8 ]
VaporPsap = ptx.StormXsection(InFname, '/all_ps_vapor', 'Vapor', 'PSAP', 
    SimCspecs, FactCspecs, 'Plots.py/FsFigVaporPsapFactors.png')
VaporPsap.Ylim = Ylim
VaporPsap.Yticks = Yticks
VaporPsap.PlotXsection()
VaporSap = ptx.StormXsection(InFname, '/all_s_vapor', 'Vapor', 'SAP', 
    SimCspecs, FactCspecs, 'Plots.py/FsFigVaporSapFactors.png')
VaporSap.Ylim = Ylim
VaporSap.Yticks = Yticks
VaporSap.PlotXsection()
