#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import FsFigContour as ffc


VelFname = 'DIAGS/ptrack_hvelocity_<SIM>.h5'
AvgFname = 'DIAGS/ptrack_avgs_<SIM>.h5'

## Ptrack Vt
#SimCspecs = [ -10, 20, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
#FactCspecs = [ -4, 4, 11 ]
#Ylim = [ 0, 10 ]
#Yticks = [ 0, 2, 4, 6, 8, 10 ]
#
#VtPsap = ffc.TrackXsection(VelFname, '/ps_v', 'Plots.py/FsFigPtrackVtPsapFactors.png', 'V_{pt}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
#VtPsap.Ylim = Ylim
#VtPsap.Yticks = Yticks
#VtPsap.CreateFig()
#
#VtSap = ffc.TrackXsection(VelFname, '/s_v', 'Plots.py/FsFigPtrackVtSapFactors.png', 'V_{pt}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
#VtSap.Yticks = Yticks
#VtSap.CreateFig()
#VtSap.CreateFig()

## Ptrack Vr
#SimCspecs = [ -5, 5, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
#FactCspecs = [ -2, 2, 11 ]
#Ylim = [ 0, 10 ]
#Yticks = [ 0, 2, 4, 6, 8, 10 ]
#
#VrPsap = ffc.TrackXsection(VelFname, '/ps_u', 'Plots.py/FsFigPtrackVrPsapFactors.png', 'V_{pr}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
#VrPsap.Ylim = Ylim
#VrPsap.Yticks = Yticks
#VrPsap.SimCmap = 'bwr'
#VrPsap.CreateFig()
#
#VrSap = ffc.TrackXsection(VelFname, '/s_u', 'Plots.py/FsFigPtrackVrSapFactors.png', 'V_{pr}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
#VrSap.Ylim = Ylim
#VrSap.Yticks = Yticks
#VrSap.SimCmap = 'bwr'
#VrSap.CreateFig()

## Theta
#SimCspecs = [ 290, 340, 11 ]
#FactCspecs = [ -3, 3, 11 ]
#Ylim = [ 0, 10 ]
#Yticks = [ 0, 2, 4, 6, 8, 10 ]
#
#ThetaPsap = ffc.TrackXsection(AvgFname, '/ps_theta', 'Plots.py/FsFigPtrackThetaPsapFactors.png', '\\theta', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
#ThetaPsap.Ylim = Ylim
#ThetaPsap.Yticks = Yticks
#ThetaPsap.CreateFig()
#
#ThetaSap = ffc.TrackXsection(AvgFname, '/s_theta', 'Plots.py/FsFigPtrackThetaSapFactors.png', '\\theta', 'SAP', SimCspecs, FactCspecs, 'ptrack')
#ThetaSap.Ylim = Ylim
#ThetaSap.Yticks = Yticks
#ThetaSap.CreateFig()

# Temp C
SimCspecs = [ -50, 30, 11 ]
FactCspecs = [ -3, 3, 11 ]
Ylim = [ 0, 10 ]
Yticks = [ 0, 2, 4, 6, 8, 10 ]

TempcPsap = ffc.TrackXsection(AvgFname, '/ps_tempc', 'Plots.py/FsFigPtrackTempcPsapFactors.png', 'T,\\ ^{\\circ}C', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
TempcPsap.Ylim = Ylim
TempcPsap.Yticks = Yticks
TempcPsap.CreateFig()

TempcSap = ffc.TrackXsection(AvgFname, '/s_tempc', 'Plots.py/FsFigPtrackTempcSapFactors.png', 'T,\\ ^{\\circ}C', 'SAP', SimCspecs, FactCspecs, 'ptrack')
TempcSap.Ylim = Ylim
TempcSap.Yticks = Yticks
TempcSap.CreateFig()

# Temp C, pressure coords
SimCspecs = [ -50, 30, 11 ]
FactCspecs = [ -3, 3, 11 ]

TempcPressPsap = ffc.TrackXsectionPress(AvgFname, '/ps_tempc_p', 'Plots.py/FsFigPtrackTempcPressPsapFactors.png', 'T,\\ ^{\\circ}C', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
TempcPressPsap.CreateFig()

TempcPressSap = ffc.TrackXsectionPress(AvgFname, '/s_tempc_p', 'Plots.py/FsFigPtrackTempcPressSapFactors.png', 'T,\\ ^{\\circ}C', 'SAP', SimCspecs, FactCspecs, 'ptrack')
TempcPressSap.CreateFig()

