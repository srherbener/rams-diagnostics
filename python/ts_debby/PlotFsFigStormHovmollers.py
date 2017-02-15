#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import FsFigContour as ffc


InFname = 'DIAGS/storm_hovs_<SIM>.h5'


## Vapor
#SimCspecs = [ 0, 20, 11 ]
#FactCspecs = [ -2, 2, 11 ]
#VaporCore = ffc.StormHovmoller(InFname, '/all_core_vapor', 'Plots.py/FsFigVaporCoreFactors.png', 'Vapor', 'CORE', SimCspecs, FactCspecs)
#VaporCore.CreateFig()
#
#VaporRband = ffc.StormHovmoller(InFname, '/all_rb_vapor', 'Plots.py/FsFigVaporRbandFactors.png', 'Vapor', 'RBAND', SimCspecs, FactCspecs)
#VaporRband.CreateFig()

## Theta
#SimCspecs = [ 290, 340, 11 ]
#FactCspecs = [ -3, 3, 11 ]
#ThetaCore = ffc.StormHovmoller(InFname, '/all_core_theta', 'Plots.py/FsFigThetaCoreFactors.png', 'Theta', 'CORE', SimCspecs, FactCspecs)
#ThetaCore.CreateFig()
#
#ThetaRband = ffc.StormHovmoller(InFname, '/all_rb_theta', 'Plots.py/FsFigThetaRbandFactors.png', 'Theta', 'RBAND', SimCspecs, FactCspecs)
#ThetaRband.CreateFig()

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

