#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import FsFigContour as ffc


InFname = 'DIAGS/storm_hovs_<SIM>.h5'


# Vapor
SimCspecs = [ 0, 20, 11 ]
FactCspecs = [ -2, 2, 11 ]
VaporCore = ffc.StormHovmoller(InFname, '/all_core_vapor', 'Plots.py/FsFigVaporCoreFactors.png', 'Vapor', 'CORE', SimCspecs, FactCspecs)
VaporCore.CreateFig()

VaporRband = ffc.StormHovmoller(InFname, '/all_rb_vapor', 'Plots.py/FsFigVaporRbandFactors.png', 'Vapor', 'RBAND', SimCspecs, FactCspecs)
VaporRband.CreateFig()

VaporLeadCore = ffc.StormHovmoller(InFname, '/lead_core_vapor', 'Plots.py/FsFigVaporLeadCoreFactors.png', 'Vapor', 'L_CORE', SimCspecs, FactCspecs)
VaporLeadCore.CreateFig()

VaporLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_vapor', 'Plots.py/FsFigVaporLeadRbandFactors.png', 'Vapor', 'L_RBAND', SimCspecs, FactCspecs)
VaporLeadRband.CreateFig()

VaporLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_vapor', 'Plots.py/FsFigVaporLeadEnvFactors.png', 'Vapor', 'L_ENV', SimCspecs, FactCspecs)
VaporLeadEnv.CreateFig()

# Theta
SimCspecs = [ 290, 340, 11 ]
FactCspecs = [ -3, 3, 11 ]
ThetaCore = ffc.StormHovmoller(InFname, '/all_core_theta', 'Plots.py/FsFigThetaCoreFactors.png', '\\theta', 'CORE', SimCspecs, FactCspecs)
ThetaCore.CreateFig()

ThetaRband = ffc.StormHovmoller(InFname, '/all_rb_theta', 'Plots.py/FsFigThetaRbandFactors.png', '\\theta', 'RBAND', SimCspecs, FactCspecs)
ThetaRband.CreateFig()

ThetaLeadCore = ffc.StormHovmoller(InFname, '/lead_core_theta', 'Plots.py/FsFigThetaLeadCoreFactors.png', '\\theta', 'L_CORE', SimCspecs, FactCspecs)
ThetaLeadCore.CreateFig()

ThetaLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_theta', 'Plots.py/FsFigThetaLeadRbandFactors.png', '\\theta', 'L_RBAND', SimCspecs, FactCspecs)
ThetaLeadRband.CreateFig()

ThetaLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_theta', 'Plots.py/FsFigThetaLeadEnvFactors.png', '\\theta', 'L_ENV', SimCspecs, FactCspecs)
ThetaLeadEnv.CreateFig()

# Theta-E
SimCspecs = [ 340, 360, 11 ]
FactCspecs = [ -8, 8, 11 ]
ThetaeCore = ffc.StormHovmoller(InFname, '/all_core_theta_e', 'Plots.py/FsFigThetaeCoreFactors.png', '\\theta_{e}', 'CORE', SimCspecs, FactCspecs)
ThetaeCore.CreateFig()

ThetaeRband = ffc.StormHovmoller(InFname, '/all_rb_theta_e', 'Plots.py/FsFigThetaeRbandFactors.png', '\\theta_{e}', 'RBAND', SimCspecs, FactCspecs)
ThetaeRband.CreateFig()

ThetaeLeadCore = ffc.StormHovmoller(InFname, '/lead_core_theta_e', 'Plots.py/FsFigThetaeLeadCoreFactors.png', '\\theta_{e}', 'L_CORE', SimCspecs, FactCspecs)
ThetaeLeadCore.CreateFig()

ThetaeLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_theta_e', 'Plots.py/FsFigThetaeLeadRbandFactors.png', '\\theta_{e}', 'L_RBAND', SimCspecs, FactCspecs)
ThetaeLeadRband.CreateFig()

ThetaeLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_theta_e', 'Plots.py/FsFigThetaeLeadEnvFactors.png', '\\theta_{e}', 'L_ENV', SimCspecs, FactCspecs)
ThetaeLeadEnv.CreateFig()

# Cloud
SimCspecs = [ 0, 0.5, 11 ]
FactCspecs = [ -0.1, 0.1, 11 ]
CloudCore = ffc.StormHovmoller(InFname, '/all_core_cloud', 'Plots.py/FsFigCloudCoreFactors.png', 'Cloud', 'CORE', SimCspecs, FactCspecs)
CloudCore.CreateFig()

CloudRband = ffc.StormHovmoller(InFname, '/all_rb_cloud', 'Plots.py/FsFigCloudRbandFactors.png', 'Cloud', 'RBAND', SimCspecs, FactCspecs)
CloudRband.CreateFig()

CloudLeadCore = ffc.StormHovmoller(InFname, '/lead_core_cloud', 'Plots.py/FsFigCloudLeadCoreFactors.png', 'Cloud', 'L_CORE', SimCspecs, FactCspecs)
CloudLeadCore.CreateFig()

CloudLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_cloud', 'Plots.py/FsFigCloudLeadRbandFactors.png', 'Cloud', 'L_RBAND', SimCspecs, FactCspecs)
CloudLeadRband.CreateFig()

CloudLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_cloud', 'Plots.py/FsFigCloudLeadEnvFactors.png', 'Cloud', 'L_ENV', SimCspecs, FactCspecs)
CloudLeadEnv.CreateFig()

# Rain
SimCspecs = [ 0, 0.5, 11 ]
FactCspecs = [ -0.1, 0.1, 11 ]
RainCore = ffc.StormHovmoller(InFname, '/all_core_rain', 'Plots.py/FsFigRainCoreFactors.png', 'Rain', 'CORE', SimCspecs, FactCspecs)
RainCore.CreateFig()

RainRband = ffc.StormHovmoller(InFname, '/all_rb_rain', 'Plots.py/FsFigRainRbandFactors.png', 'Rain', 'RBAND', SimCspecs, FactCspecs)
RainRband.CreateFig()

RainLeadCore = ffc.StormHovmoller(InFname, '/lead_core_rain', 'Plots.py/FsFigRainLeadCoreFactors.png', 'Rain', 'L_CORE', SimCspecs, FactCspecs)
RainLeadCore.CreateFig()

RainLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_rain', 'Plots.py/FsFigRainLeadRbandFactors.png', 'Rain', 'L_RBAND', SimCspecs, FactCspecs)
RainLeadRband.CreateFig()

RainLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_rain', 'Plots.py/FsFigRainLeadEnvFactors.png', 'Rain', 'L_ENV', SimCspecs, FactCspecs)
RainLeadEnv.CreateFig()

# Cloud Evap
SimCspecs = [ -2.0, 0.0, 11 ]
FactCspecs = [ -0.5, 0.5, 11 ]
CloudEvapCore = ffc.StormHovmoller(InFname, '/all_core_cloud_evap', 'Plots.py/FsFigCloudEvapCoreFactors.png', 'Cloud Evap', 'CORE', SimCspecs, FactCspecs)
CloudEvapCore.CreateFig()

CloudEvapRband = ffc.StormHovmoller(InFname, '/all_rb_cloud_evap', 'Plots.py/FsFigCloudEvapRbandFactors.png', 'Cloud Evap', 'RBAND', SimCspecs, FactCspecs)
CloudEvapRband.CreateFig()

CloudEvapLeadCore = ffc.StormHovmoller(InFname, '/lead_core_cloud_evap', 'Plots.py/FsFigCloudEvapLeadCoreFactors.png', 'Cloud Evap', 'L_CORE', SimCspecs, FactCspecs)
CloudEvapLeadCore.CreateFig()

CloudEvapLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_cloud_evap', 'Plots.py/FsFigCloudEvapLeadRbandFactors.png', 'Cloud Evap', 'L_RBAND', SimCspecs, FactCspecs)
CloudEvapLeadRband.CreateFig()

CloudEvapLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_cloud_evap', 'Plots.py/FsFigCloudEvapLeadEnvFactors.png', 'Cloud Evap', 'L_ENV', SimCspecs, FactCspecs)
CloudEvapLeadEnv.CreateFig()

# Cloud Cond
SimCspecs = [ 0.0, 2.0, 11 ]
FactCspecs = [ -0.5, 0.5, 11 ]
CloudCondCore = ffc.StormHovmoller(InFname, '/all_core_cloud_cond', 'Plots.py/FsFigCloudCondCoreFactors.png', 'Cloud Cond', 'CORE', SimCspecs, FactCspecs)
CloudCondCore.CreateFig()

CloudCondRband = ffc.StormHovmoller(InFname, '/all_rb_cloud_cond', 'Plots.py/FsFigCloudCondRbandFactors.png', 'Cloud Cond', 'RBAND', SimCspecs, FactCspecs)
CloudCondRband.CreateFig()

CloudCondLeadCore = ffc.StormHovmoller(InFname, '/lead_core_cloud_cond', 'Plots.py/FsFigCloudCondLeadCoreFactors.png', 'Cloud Cond', 'L_CORE', SimCspecs, FactCspecs)
CloudCondLeadCore.CreateFig()

CloudCondLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_cloud_cond', 'Plots.py/FsFigCloudCondLeadRbandFactors.png', 'Cloud Cond', 'L_RBAND', SimCspecs, FactCspecs)
CloudCondLeadRband.CreateFig()

CloudCondLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_cloud_cond', 'Plots.py/FsFigCloudCondLeadEnvFactors.png', 'Cloud Cond', 'L_ENV', SimCspecs, FactCspecs)
CloudCondLeadEnv.CreateFig()

# Rain Evap
SimCspecs = [ -1.0, 0.0, 11 ]
FactCspecs = [ -0.5, 0.5, 11 ]
RainEvapCore = ffc.StormHovmoller(InFname, '/all_core_rain_evap', 'Plots.py/FsFigRainEvapCoreFactors.png', 'Rain Evap', 'CORE', SimCspecs, FactCspecs)
RainEvapCore.CreateFig()
 
RainEvapRband = ffc.StormHovmoller(InFname, '/all_rb_rain_evap', 'Plots.py/FsFigRainEvapRbandFactors.png', 'Rain Evap', 'RBAND', SimCspecs, FactCspecs)
RainEvapRband.CreateFig()

RainEvapLeadCore = ffc.StormHovmoller(InFname, '/lead_core_rain_evap', 'Plots.py/FsFigRainEvapLeadCoreFactors.png', 'Rain Evap', 'L_CORE', SimCspecs, FactCspecs)
RainEvapLeadCore.CreateFig()

RainEvapLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_rain_evap', 'Plots.py/FsFigRainEvapLeadRbandFactors.png', 'Rain Evap', 'L_RBAND', SimCspecs, FactCspecs)
RainEvapLeadRband.CreateFig()

RainEvapLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_rain_evap', 'Plots.py/FsFigRainEvapLeadEnvFactors.png', 'Rain Evap', 'L_ENV', SimCspecs, FactCspecs)
RainEvapLeadEnv.CreateFig()

# Rain Cond
SimCspecs = [ 0.0, 0.5, 11 ]
FactCspecs = [ -0.1, 0.1, 11 ]
RainCondCore = ffc.StormHovmoller(InFname, '/all_core_rain_cond', 'Plots.py/FsFigRainCondCoreFactors.png', 'Rain Cond', 'CORE', SimCspecs, FactCspecs)
RainCondCore.CreateFig()
 
RainCondRband = ffc.StormHovmoller(InFname, '/all_rb_rain_cond', 'Plots.py/FsFigRainCondRbandFactors.png', 'Rain Cond', 'RBAND', SimCspecs, FactCspecs)
RainCondRband.CreateFig()

RainCondLeadCore = ffc.StormHovmoller(InFname, '/lead_core_rain_cond', 'Plots.py/FsFigRainCondLeadCoreFactors.png', 'Rain Cond', 'L_CORE', SimCspecs, FactCspecs)
RainCondLeadCore.CreateFig()

RainCondLeadRband = ffc.StormHovmoller(InFname, '/lead_rb_rain_cond', 'Plots.py/FsFigRainCondLeadRbandFactors.png', 'Rain Cond', 'L_RBAND', SimCspecs, FactCspecs)
RainCondLeadRband.CreateFig()

RainCondLeadEnv = ffc.StormHovmoller(InFname, '/lead_env_rain_cond', 'Plots.py/FsFigRainCondLeadEnvFactors.png', 'Rain Cond', 'L_ENV', SimCspecs, FactCspecs)
RainCondLeadEnv.CreateFig()

# IceDep
SimCspecs = [ 0, 1.5, 11 ]
FactCspecs = [ -1, 1, 11 ]
Ylim = [ 0, 15 ]
Yticks = [ 0, 3, 6, 9, 12, 15 ]

IceDepCore = ffc.StormHovmoller(InFname, '/all_core_ice_dep', 'Plots.py/FsFigIceDepCoreFactors.png', 'Ice Dep', 'CORE', SimCspecs, FactCspecs)
IceDepCore.Ylim = Ylim
IceDepCore.Yticks = Yticks
IceDepCore.CreateFig()

IceDepRband = ffc.StormHovmoller(InFname, '/all_rb_ice_dep', 'Plots.py/FsFigIceDepRbandFactors.png', 'Ice Dep', 'RBAND', SimCspecs, FactCspecs)
IceDepRband.Ylim = Ylim
IceDepRband.Yticks = Yticks
IceDepRband.CreateFig()

# Rain to ice
SimCspecs = [ 0, 0.3, 11 ]
FactCspecs = [ -0.15, 0.15, 11 ]
Ylim = [ 0, 15 ]
Yticks = [ 0, 3, 6, 9, 12, 15 ]

Rain2IceCore = ffc.StormHovmoller(InFname, '/all_core_rain2ice', 'Plots.py/FsFigRain2IceCoreFactors.png', 'Rain2Ice', 'CORE', SimCspecs, FactCspecs)
Rain2IceCore.CreateFig()

Rain2IceRband = ffc.StormHovmoller(InFname, '/all_rb_rain2ice', 'Plots.py/FsFigRain2IceRbandFactors.png', 'Rain2Ice', 'RBAND', SimCspecs, FactCspecs)
Rain2IceRband.CreateFig()

# Updrafts
SimCspecs = [ 0, 0.3, 11 ]
FactCspecs = [ -0.1, 0.1, 11 ]
Ylim = [ 0, 15 ]
Yticks = [ 0, 3, 6, 9, 12, 15 ]

UpdraftCore = ffc.StormHovmoller(InFname, '/all_core_updraft', 'Plots.py/FsFigUpdraftCoreFactors.png', 'w_{up}', 'CORE', SimCspecs, FactCspecs)
UpdraftCore.Ylim = Ylim
UpdraftCore.Yticks = Yticks
UpdraftCore.CreateFig()

SimCspecs = [ 0, 0.2, 11 ]
FactCspecs = [ -0.05, 0.05, 11 ]
UpdraftRband = ffc.StormHovmoller(InFname, '/all_rb_updraft', 'Plots.py/FsFigUpdraftRbandFactors.png', 'w_{up}', 'RBAND', SimCspecs, FactCspecs)
UpdraftRband.Ylim = Ylim
UpdraftRband.Yticks = Yticks
UpdraftRband.CreateFig()

# Downdrafts
SimCspecs = [ -0.2, 0.0, 11 ]
FactCspecs = [ -0.05, 0.05, 11 ]
Ylim = [ 0, 15 ]
Yticks = [ 0, 3, 6, 9, 12, 15 ]

DowndraftCore = ffc.StormHovmoller(InFname, '/all_core_dndraft', 'Plots.py/FsFigDowndraftCoreFactors.png', 'w_{dn}', 'CORE', SimCspecs, FactCspecs)
DowndraftCore.Ylim = Ylim
DowndraftCore.Yticks = Yticks
DowndraftCore.CreateFig()

SimCspecs = [ -0.15, 0.0, 11 ]
FactCspecs = [ -0.02, 0.02, 11 ]
DowndraftRband = ffc.StormHovmoller(InFname, '/all_rb_dndraft', 'Plots.py/FsFigDowndraftRbandFactors.png', 'w_{dn}', 'RBAND', SimCspecs, FactCspecs)
DowndraftRband.Ylim = Ylim
DowndraftRband.Yticks = Yticks
DowndraftRband.CreateFig()

