#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import FsFigContour as ffc


InFname = 'DIAGS/storm_xsections_<SIM>.h5'

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
## Vr
#SimCspecs = [ -5, 5, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
#FactCspecs = [ -2, 2, 21 ]
#Ylim = [ 0, 10 ]
#Yticks = [ 0, 2, 4, 6, 8, 10 ]
#VrPsap = ptx.StormXsection(InFname, '/all_ps_speed_r', 'V_r', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigVrPsapFactors.png')
#VrPsap.Ylim = Ylim
#VrPsap.Yticks = Yticks
#VrPsap.Cmap = 'bwr'
#VrPsap.PlotXsection()
#VrSap = ptx.StormXsection(InFname, '/all_s_speed_r', 'V_r', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigVrSapFactors.png')
#VrSap.Ylim = Ylim
#VrSap.Yticks = Yticks
#VrSap.Cmap = 'bwr'
#VrSap.PlotXsection()
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
#SimCspecs = [ 340, 360, 11 ]
#FactCspecs = [ -8, 8, 11 ]
#Ylim = [ 0, 8 ]
#Yticks =  [ 0, 2, 4, 6, 8 ]
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
#
## Theta
#SimCspecs = [ 290, 340, 11 ]
#FactCspecs = [ -3, 3, 11 ]
#Ylim = [ 0, 8 ]
#Yticks =  [ 0, 2, 4, 6, 8 ]
#ThetaPsap = ptx.StormXsection(InFname, '/all_ps_theta', 'Theta', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigThetaPsapFactors.png')
#ThetaPsap.Ylim = Ylim
#ThetaPsap.Yticks = Yticks
#ThetaPsap.PlotXsection()
#ThetaSap = ptx.StormXsection(InFname, '/all_s_theta', 'Theta', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigThetaSapFactors.png')
#ThetaSap.Ylim = Ylim
#ThetaSap.Yticks = Yticks
#ThetaSap.PlotXsection()
#
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
#
## Vapor
#SimCspecs = [ 0, 20, 11 ]
#FactCspecs = [ -2, 2, 11 ]
#Ylim = [ 0, 8 ]
#Yticks =  [ 0, 2, 4, 6, 8 ]
#VaporPsap = ptx.StormXsection(InFname, '/all_ps_vapor', 'Vapor', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigVaporPsapFactors.png')
#VaporPsap.Ylim = Ylim
#VaporPsap.Yticks = Yticks
#VaporPsap.PlotXsection()
#VaporSap = ptx.StormXsection(InFname, '/all_s_vapor', 'Vapor', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigVaporSapFactors.png')
#VaporSap.Ylim = Ylim
#VaporSap.Yticks = Yticks
#VaporSap.PlotXsection()
#
## Vapor (in lead region)
#SimCspecs = [ 0, 20, 11 ]
#FactCspecs = [ -2, 2, 11 ]
#Ylim = [ 0, 8 ]
#Yticks =  [ 0, 2, 4, 6, 8 ]
#VaporLeadPsap = ptx.StormXsection(InFname, '/lead_ps_vapor', 'VaporLead', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigVaporLeadPsapFactors.png')
#VaporLeadPsap.Ylim = Ylim
#VaporLeadPsap.Yticks = Yticks
#VaporLeadPsap.PlotXsection()
#VaporLeadSap = ptx.StormXsection(InFname, '/lead_s_vapor', 'VaporLead', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigVaporLeadSapFactors.png')
#VaporLeadSap.Ylim = Ylim
#VaporLeadSap.Yticks = Yticks
#VaporLeadSap.PlotXsection()
#
## Cloud
#SimCspecs = [ 0, 0.5, 11 ]
#FactCspecs = [ -0.1, 0.1, 11 ]
#Ylim = [ 0, 8 ]
#Yticks =  [ 0, 2, 4, 6, 8 ]
#CloudPsap = ptx.StormXsection(InFname, '/all_ps_cloud_mass', 'Cloud', 'PSAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigCloudPsapFactors.png')
#CloudPsap.Ylim = Ylim
#CloudPsap.Yticks = Yticks
#CloudPsap.PlotXsection()
#CloudSap = ptx.StormXsection(InFname, '/all_s_cloud_mass', 'Cloud', 'SAP', 
#    SimCspecs, FactCspecs, 'Plots.py/FsFigCloudSapFactors.png')
#CloudSap.Ylim = Ylim
#CloudSap.Yticks = Yticks
#CloudSap.PlotXsection()

# Rain
SimCspecs = [    0, 0.5, 11 ]
FacCspecs = [ -0.1, 0.1, 11 ]
RainPsap = ffc.StormXsection(InFname, '/all_ps_rain_mass', 'Plots.py/FsFigRainPsapFactors.png', 'Rain', 'PSAP', SimCspecs, FacCspecs)
RainPsap.CreateFig()

RainSap = ffc.StormXsection(InFname, '/all_s_rain_mass', 'Plots.py/FsFigRainSapFactors.png', 'Rain', 'SAP', SimCspecs, FacCspecs)
RainSap.CreateFig()
#
