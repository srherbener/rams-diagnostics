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
#Ylim = [ 0, 10 ]
#Yticks = [ 0, 2, 4, 6, 8, 10 ]
#
#VtPsap = ffc.StormXsection(InFname, '/all_ps_speed_t', 'Plots.py/FsFigVtPsapFactors.png', 'V_t', 'PSAP', SimCspecs, FactCspecs)
#VtPsap.Ylim = Ylim
#VtPsap.Yticks = Yticks
#VtPsap.CreateFig()
#
#VtSap = ffc.StormXsection(InFname, '/all_s_speed_t', 'Plots.py/FsFigVtSapFactors.png', 'V_t', 'SAP', SimCspecs, FactCspecs)
#VtSap.Yticks = Yticks
#VtSap.CreateFig()
#
## Vr
#SimCspecs = [ -5, 5, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
#FactCspecs = [ -2, 2, 21 ]
#Ylim = [ 0, 10 ]
#Yticks = [ 0, 2, 4, 6, 8, 10 ]
#
#VrPsap = ffc.StormXsection(InFname, '/all_ps_speed_r', 'Plots.py/FsFigVrPsapFactors.png', 'V_r', 'PSAP', SimCspecs, FactCspecs)
#VrPsap.Ylim = Ylim
#VrPsap.Yticks = Yticks
#VrPsap.SimCmap = 'bwr'
#VrPsap.CreateFig()
#
#VrSap = ffc.StormXsection(InFname, '/all_s_speed_r', 'Plots.py/FsFigVrSapFactors.png', 'V_r', 'SAP', SimCspecs, FactCspecs)
#VrSap.Ylim = Ylim
#VrSap.Yticks = Yticks
#VrSap.SimCmap = 'bwr'
#VrSap.CreateFig()
#
## Theta-E
#SimCspecs = [ 340, 360, 11 ]
#FactCspecs = [ -8, 8, 11 ]
#ThetaePsap = ffc.StormXsection(InFname, '/all_ps_theta_e', 'Plots.py/FsFigThetaePsapFactors.png', '\\theta_{e}', 'PSAP', SimCspecs, FactCspecs)
#ThetaePsap.CreateFig()
#ThetaeSap = ffc.StormXsection(InFname, '/all_s_theta_e', 'Plots.py/FsFigThetaeSapFactors.png', '\\theta_{e}', 'SAP', SimCspecs, FactCspecs)
#ThetaeSap.CreateFig()
#
## Vapor
#SimCspecs = [ 0, 20, 11 ]
#FactCspecs = [ -2, 2, 11 ]
#VaporPsap = ffc.StormXsection(InFname, '/all_ps_vapor', 'Plots.py/FsFigVaporPsapFactors.png', 'Vapor', 'PSAP', SimCspecs, FactCspecs)
#VaporPsap.CreateFig()
#VaporSap = ffc.StormXsection(InFname, '/all_s_vapor', 'Plots.py/FsFigVaporSapFactors.png', 'Vapor', 'SAP', SimCspecs, FactCspecs)
#VaporSap.CreateFig()
#
## Cloud
#SimCspecs = [ 0, 0.5, 11 ]
#FactCspecs = [ -0.1, 0.1, 11 ]
#CloudPsap = ffc.StormXsection(InFname,  '/all_ps_cloud_mass', 'Plots.py/FsFigCloudPsapFactors.png', 'Cloud', 'PSAP', SimCspecs, FactCspecs)
#CloudPsap.CreateFig()
#CloudSap = ffc.StormXsection(InFname, '/all_s_cloud_mass', 'Plots.py/FsFigCloudSapFactors.png', 'Cloud', 'SAP', SimCspecs, FactCspecs)
#CloudSap.CreateFig()
#
## Rain
#SimCspecs = [    0, 0.5, 11 ]
#FacCspecs = [ -0.1, 0.1, 11 ]
#RainPsap = ffc.StormXsection(InFname, '/all_ps_rain_mass', 'Plots.py/FsFigRainPsapFactors.png', 'Rain', 'PSAP', SimCspecs, FacCspecs)
#RainPsap.CreateFig()
#RainSap = ffc.StormXsection(InFname, '/all_s_rain_mass', 'Plots.py/FsFigRainSapFactors.png', 'Rain', 'SAP', SimCspecs, FacCspecs)
#RainSap.CreateFig()
#
## Ice deposition
#SimCspecs = [ 0, 1.5, 16 ]
#FactCspecs = [ -1, 1, 11 ]
#Ylim = [ 0, 15 ]
#Yticks = [ 0, 3, 6, 9, 12, 15 ]
#IceDepPsap = ffc.StormXsection(InFname, '/all_ps_ice_dep', 'Plots.py/FsFigIceDepPsapFactors.png', 'Ice Dep', 'PSAP', SimCspecs, FactCspecs)
#IceDepPsap.Ylim = Ylim
#IceDepPsap.Yticks = Yticks
#IceDepPsap.CreateFig()
#IceDepSap = ffc.StormXsection(InFname, '/all_s_ice_dep', 'Plots.py/FsFigIceDepSapFactors.png', 'Ice Dep', 'SAP', SimCspecs, FactCspecs)
#IceDepSap.Ylim = Ylim
#IceDepSap.Yticks = Yticks
#IceDepSap.CreateFig()

##### rime
####SimCspecs = [ 0, 1, 11 ]
####FactCspecs = [ -1, 1, 11 ]
####Ylim = [ 0, 12 ]
####Yticks =  [ 0, 3, 6, 9, 12 ]
####RimePsap = ffc.StormXsection(InFname, '/all_ps_cloud_rime', 'Plots.py/FsFigRimePsapFactors.png', 'Rime', 'PSAP', SimCspecs, FactCspecs)
####RimePsap.Ylim = Ylim
####RimePsap.Yticks = Yticks
####RimePsap.CreateFig()
####RimeSap = ffc.StormXsection(InFname, '/all_s_cloud_rime', 'Plots.py/FsFigRimeSapFactors.png', 'Rime', 'SAP', SimCspecs, FactCspecs)
####RimeSap.Ylim = Ylim
####RimeSap.Yticks = Yticks
####RimeSap.CreateFig()

##### rain to ice
####SimCspecs = [ 0, 0.3, 11 ]
####FactCspecs = [ -0.15, 0.15, 11 ]
####Ylim = [ 0, 12 ]
####Yticks =  [ 0, 3, 6, 9, 12 ]
####R2icePsap = ffc.StormXsection(InFname, '/all_ps_rain2ice', 'Plots.py/FsFigR2icePsapFactors.png', 'Rain2Ice', 'PSAP', SimCspecs, FactCspecs)
####R2icePsap.Ylim = Ylim
####R2icePsap.Yticks = Yticks
####R2icePsap.CreateFig()
####R2iceSap = ffc.StormXsection(InFname, '/all_s_rain2ice', 'Plots.py/FsFigR2iceSapFactors.png', 'Rain2Ice', 'SAP', SimCspecs, FactCspecs)
####R2iceSap.Ylim = Ylim
####R2iceSap.Yticks = Yticks
####R2iceSap.CreateFig()

##### Theta
####SimCspecs = [ 290, 340, 11 ]
####FactCspecs = [ -3, 3, 11 ]
####ThetaPsap = ffc.StormXsection(InFname, '/all_ps_theta', 'Plots.py/FsFigThetaPsapFactors.png', '\\theta', 'PSAP', SimCspecs, FactCspecs)
####ThetaPsap.CreateFig()
####ThetaSap = ffc.StormXsection(InFname, '/all_s_theta', 'Plots.py/FsFigThetaSapFactors.png', '\\theta', 'SAP', SimCspecs, FactCspecs)
####ThetaSap.CreateFig()

##### Cloud Evap
####SimCspecs = [ -2, 0, 11 ]
####FactCspecs = [ -0.5, 0.5, 11 ]
####CloudEvapPsap = ffc.StormXsection(InFname, '/all_ps_cloud_evap', 'Plots.py/FsFigCloudEvapPsapFactors.png', 'Cloud Evap', 'PSAP', SimCspecs, FactCspecs)
####CloudEvapPsap.CreateFig()
####CloudEvapSap = ffc.StormXsection(InFname, '/all_s_cloud_evap', 'Plots.py/FsFigCloudEvapSapFactors.png', 'Cloud Evap', 'SAP', SimCspecs, FactCspecs)
####CloudEvapSap.CreateFig()

##### Rain Evap
####SimCspecs = [ -1, 0, 11 ]
####FactCspecs = [ -0.5, 0.5, 11 ]
####RainEvapPsap = ffc.StormXsection(InFname, '/all_ps_rain_evap', 'Plots.py/FsFigRainEvapPsapFactors.png', 'Rain Evap', 'PSAP', SimCspecs, FactCspecs)
####RainEvapPsap.CreateFig()
####RainEvapSap = ffc.StormXsection(InFname, '/all_s_rain_evap', 'Plots.py/FsFigRainEvapSapFactors.png', 'Rain Evap', 'SAP', SimCspecs, FactCspecs)
####RainEvapSap.CreateFig()

##### Cloud Cond
####SimCspecs = [ 0, 2, 11 ]
####FactCspecs = [ -0.5, 0.5, 11 ]
####CloudCondPsap = ffc.StormXsection(InFname, '/all_ps_cloud_cond', 'Plots.py/FsFigCloudCondPsap.png', 'Cloud Cond', 'PSAP', SimCspecs, FactCspecs)
####CloudCondPsap.CreateFig()
####CloudCondSap = ffc.StormXsection(InFname, '/all_s_cloud_cond', 'Plots.py/FsFigCloudCondSap.png', 'Cloud Cond', 'SAP', SimCspecs, FactCspecs)
####CloudCondSap.CreateFig()

##### Rain Cond
####SimCspecs = [ 0, 0.5, 11 ]
####FactCspecs = [ -0.1, 0.1, 11 ]
####RainCondPsap = ffc.StormXsection(InFname, '/all_ps_rain_cond', 'Plots.py/FsFigRainCondPsap.png', 'Rain Cond', 'PSAP', SimCspecs, FactCspecs)
####RainCondPsap.CreateFig()
####RainCondSap = ffc.StormXsection(InFname, '/all_s_rain_cond', 'Plots.py/FsFigRainCondSap.png', 'Rain Cond', 'SAP', SimCspecs, FactCspecs)
####RainCondSap.CreateFig()

##### Vapor (in lead region)
####SimCspecs = [ 0, 20, 11 ]
####FactCspecs = [ -2, 2, 11 ]
####VaporLeadPsap = ffc.StormXsection(InFname, '/lead_ps_vapor', 'Plots.py/FsFigVaporLeadPsapFactors.png', 'VaporLead', 'PSAP', SimCspecs, FactCspecs)
####VaporLeadPsap.CreateFig()
####VaporLeadSap = ffc.StormXsection(InFname, '/lead_s_vapor', 'Plots.py/FsFigVaporLeadSapFactors.png', 'VaporLead', 'SAP', SimCspecs, FactCspecs)
####VaporLeadSap.CreateFig()

############################################################################################
# Pressure coords
############################################################################################

## Max Vt
#SimCspecs = [ 0, 20, 21 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
#FactCspecs = [ -4, 4, 21 ]
#
#VtPsap = ffc.StormXsectionPress(InFname, '/all_p_ps_speed_t', 'Plots.py/FsFigVtPressPsapFactors.png', 'V_t', 'PSAP', SimCspecs, FactCspecs)
#VtPsap.CreateFig()
#
#VtSap = ffc.StormXsectionPress(InFname, '/all_p_s_speed_t', 'Plots.py/FsFigVtPressSapFactors.png', 'V_t', 'SAP', SimCspecs, FactCspecs)
#VtSap.CreateFig()

## Vr
#SimCspecs = [ -5, 5, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
#FactCspecs = [ -2, 2, 21 ]
#
#VrPsap = ffc.StormXsectionPress(InFname, '/all_p_ps_speed_r', 'Plots.py/FsFigVrPressPsapFactors.png', 'V_r', 'PSAP', SimCspecs, FactCspecs)
#VrPsap.SimCmap = 'bwr'
#VrPsap.CreateFig()
#
#VrSap = ffc.StormXsectionPress(InFname, '/all_p_s_speed_r', 'Plots.py/FsFigVrPressSapFactors.png', 'V_r', 'SAP', SimCspecs, FactCspecs)
#VrSap.SimCmap = 'bwr'
#VrSap.CreateFig()

# Updrafts
#SimCspecs = [ 0, 0.3, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
#FactCspecs = [ -0.1, 0.1, 21 ]
#
#UpPsap = ffc.StormXsectionPress(InFname, '/all_p_ps_updraft', 'Plots.py/FsFigUpPressPsapFactors.png', 'W_{up}', 'PSAP', SimCspecs, FactCspecs)
#UpPsap.CreateFig()
#
#UpSap = ffc.StormXsectionPress(InFname, '/all_p_s_updraft', 'Plots.py/FsFigUpPressSapFactors.png', 'W_{up}', 'SAP', SimCspecs, FactCspecs)
#UpSap.CreateFig()

## Theta-E
#SimCspecs = [ 330, 350, 11 ]
#FactCspecs = [ -2, 2, 11 ]
#ThetaePsap = ffc.StormXsectionPress(InFname, '/all_p_ps_theta_e', 'Plots.py/FsFigThetaePressPsapFactors.png', '\\theta_{e}', 'PSAP', SimCspecs, FactCspecs)
#ThetaePsap.CreateFig()
#ThetaeSap = ffc.StormXsectionPress(InFname, '/all_p_s_theta_e', 'Plots.py/FsFigThetaePressSapFactors.png', '\\theta_{e}', 'SAP', SimCspecs, FactCspecs)
#ThetaeSap.CreateFig()

## Theta-E, vortex removed
#SimCspecs = [ 330, 350, 11 ]
#FactCspecs = [ -2, 2, 11 ]
#ThetaePsap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_ps_theta_e_lite', 'Plots.py/FsFigThetaeNvPressPsapFactors.png', '\\theta_{e}', 'PSAP', SimCspecs, FactCspecs)
#ThetaePsap.CreateFig()
#ThetaeSap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_s_theta_e_lite', 'Plots.py/FsFigThetaeNvPressSapFactors.png', '\\theta_{e}', 'SAP', SimCspecs, FactCspecs)
#ThetaeSap.CreateFig()

## Vapor
#SimCspecs = [ 0, 20, 11 ]
#FactCspecs = [ -1, 1, 11 ]
#VaporPsap = ffc.StormXsectionPressLite(InFname, '/all_p_ps_vapor_lite', 'Plots.py/FsFigVaporPressPsapFactors.png', 'R_{v}', 'PSAP', SimCspecs, FactCspecs)
#VaporPsap.CreateFig()
#VaporSap = ffc.StormXsectionPressLite(InFname, '/all_p_s_vapor_lite', 'Plots.py/FsFigVaporPressSapFactors.png', 'R_{v}', 'SAP', SimCspecs, FactCspecs)
#VaporSap.CreateFig()

## Vapor, vortex removed
#SimCspecs = [ 0, 20, 11 ]
#FactCspecs = [ -1, 1, 11 ]
#VaporPsap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_ps_vapor_lite', 'Plots.py/FsFigVaporNvPressPsapFactors.png', 'R_{v}', 'PSAP', SimCspecs, FactCspecs)
#VaporPsap.CreateFig()
#VaporSap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_s_vapor_lite', 'Plots.py/FsFigVaporNvPressSapFactors.png', 'R_{v}', 'SAP', SimCspecs, FactCspecs)
#VaporSap.CreateFig()

## Theta
#SimCspecs = [ 300, 350, 11 ]
#FactCspecs = [ -2, 2, 11 ]
#ThetaPsap = ffc.StormXsectionPressLite(InFname, '/all_p_ps_theta_lite', 'Plots.py/FsFigThetaPressPsapFactors.png', '\\theta', 'PSAP', SimCspecs, FactCspecs)
#ThetaPsap.CreateFig()
#ThetaSap = ffc.StormXsectionPressLite(InFname, '/all_p_s_theta_lite', 'Plots.py/FsFigThetaPressSapFactors.png', '\\theta', 'SAP', SimCspecs, FactCspecs)
#ThetaSap.CreateFig()

## Theta, vortex removed
#SimCspecs = [ 300, 350, 11 ]
#FactCspecs = [ -2, 2, 11 ]
#ThetaPsap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_ps_theta_lite', 'Plots.py/FsFigThetaNvPressPsapFactors.png', '\\theta', 'PSAP', SimCspecs, FactCspecs)
#ThetaPsap.CreateFig()
#ThetaSap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_s_theta_lite', 'Plots.py/FsFigThetaNvPressSapFactors.png', '\\theta', 'SAP', SimCspecs, FactCspecs)
#ThetaSap.CreateFig()

# LiqEvap
SimCspecs = [ -0.5, 0, 11 ]
FactCspecs = [ -0.1, 0.1, 11 ]
LiqEvapPsap = ffc.StormXsectionPressLite(InFname, '/all_p_ps_liq_evap_lite', 'Plots.py/FsFigLiqEvapPressPsapFactors.png', 'Evap', 'PSAP', SimCspecs, FactCspecs)
LiqEvapPsap.CreateFig()
LiqEvapSap = ffc.StormXsectionPressLite(InFname, '/all_p_s_liq_evap_lite', 'Plots.py/FsFigLiqEvapPressSapFactors.png', 'Evap', 'SAP', SimCspecs, FactCspecs)
LiqEvapSap.CreateFig()

# LiqEvap, vortex removed
SimCspecs = [ -0.5, 0, 11 ]
FactCspecs = [ -0.1, 0.1, 11 ]
LiqEvapPsap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_ps_liq_evap_lite', 'Plots.py/FsFigLiqEvapNvPressPsapFactors.png', 'Evap', 'PSAP', SimCspecs, FactCspecs)
LiqEvapPsap.CreateFig()
LiqEvapSap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_s_liq_evap_lite', 'Plots.py/FsFigLiqEvapNvPressSapFactors.png', 'Evap', 'SAP', SimCspecs, FactCspecs)
LiqEvapSap.CreateFig()

# LiqCond
SimCspecs = [ 0, 0.5, 11 ]
FactCspecs = [ -0.05, 0.05, 11 ]
LiqCondPsap = ffc.StormXsectionPressLite(InFname, '/all_p_ps_liq_cond_lite', 'Plots.py/FsFigLiqCondPressPsapFactors.png', 'Cond', 'PSAP', SimCspecs, FactCspecs)
LiqCondPsap.CreateFig()
LiqCondSap = ffc.StormXsectionPressLite(InFname, '/all_p_s_liq_cond_lite', 'Plots.py/FsFigLiqCondPressSapFactors.png', 'Cond', 'SAP', SimCspecs, FactCspecs)
LiqCondSap.CreateFig()

# LiqCond, vortex removed
SimCspecs = [ 0, 0.5, 11 ]
FactCspecs = [ -0.05, 0.05, 11 ]
LiqCondPsap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_ps_liq_cond_lite', 'Plots.py/FsFigLiqCondNvPressPsapFactors.png', 'Cond', 'PSAP', SimCspecs, FactCspecs)
LiqCondPsap.CreateFig()
LiqCondSap = ffc.StormXsectionPressLite(InFname, '/all_p_nv_s_liq_cond_lite', 'Plots.py/FsFigLiqCondNvPressSapFactors.png', 'Cond', 'SAP', SimCspecs, FactCspecs)
LiqCondSap.CreateFig()

