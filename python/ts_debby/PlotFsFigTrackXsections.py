#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import FsFigContour as ffc


VelFname = 'DIAGS/ptrack_hvel_<SIM>.h5'
AvgFname = 'DIAGS/ptrack_avgs_<SIM>.h5'

# Ptrack Vt
SimCspecs = [ -10, 20, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -4, 4, 11 ]

VtPsap = ffc.TrackXsection(VelFname, '/ps_v', 'Plots.py/FsFigPtrackVtPsapFactors.png', 'V_{pt}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
VtPsap.CreateFig()

VtSap = ffc.TrackXsection(VelFname, '/s_v', 'Plots.py/FsFigPtrackVtSapFactors.png', 'V_{pt}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
VtSap.CreateFig()

# Ptrack Vr
SimCspecs = [ -5, 5, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -2, 2, 11 ]

VrPsap = ffc.TrackXsection(VelFname, '/ps_u', 'Plots.py/FsFigPtrackVrPsapFactors.png', 'V_{pr}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
VrPsap.SimCmap = 'bwr'
VrPsap.CreateFig()

VrSap = ffc.TrackXsection(VelFname, '/s_u', 'Plots.py/FsFigPtrackVrSapFactors.png', 'V_{pr}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
VrSap.SimCmap = 'bwr'
VrSap.CreateFig()

# Ptrack Vt, pressure coords
SimCspecs = [ -10, 20, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -4, 4, 11 ]

VtPressPsap = ffc.TrackXsectionPress(VelFname, '/ps_v_p', 'Plots.py/FsFigPtrackVtPressPsapFactors.png', 'V_{pt}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
VtPressPsap.CreateFig()

VtPressSap = ffc.TrackXsectionPress(VelFname, '/s_v_p', 'Plots.py/FsFigPtrackVtPressSapFactors.png', 'V_{pt}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
VtPressSap.CreateFig()

# Ptrack Vr, pressure coords
SimCspecs = [ -5, 5, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -2, 2, 11 ]

VrPressPsap = ffc.TrackXsectionPress(VelFname, '/ps_u_p', 'Plots.py/FsFigPtrackVrPressPsapFactors.png', 'V_{pr}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
VrPressPsap.SimCmap = 'bwr'
VrPressPsap.CreateFig()

VrPressSap = ffc.TrackXsectionPress(VelFname, '/s_u_p', 'Plots.py/FsFigPtrackVrPressSapFactors.png', 'V_{pr}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
VrPressSap.SimCmap = 'bwr'
VrPressSap.CreateFig()


# Ptrack Vt, vortex removed
SimCspecs = [ -10, 20, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -4, 4, 11 ]
Xname = '/xl_coords'

VtNvPsap = ffc.TrackXsection(VelFname, '/ps_v_nv_lite', 'Plots.py/FsFigPtrackVtNvPsapFactors.png', 'V_{ptnv}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
VtNvPsap.XcoordName = Xname
VtNvPsap.CreateFig()

VtNvSap = ffc.TrackXsection(VelFname, '/s_v_nv_lite', 'Plots.py/FsFigPtrackVtNvSapFactors.png', 'V_{ptnv}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
VtNvSap.XcoordName = Xname
VtNvSap.CreateFig()

# Ptrack Vr, vortex removed
SimCspecs = [ -5, 5, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -2, 2, 11 ]

VrNvPsap = ffc.TrackXsection(VelFname, '/ps_u_nv_lite', 'Plots.py/FsFigPtrackVrNvPsapFactors.png', 'V_{prnv}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
VrNvPsap.XcoordName = Xname
VrNvPsap.SimCmap = 'bwr'
VrNvPsap.CreateFig()

VrNvSap = ffc.TrackXsection(VelFname, '/s_u_nv_lite', 'Plots.py/FsFigPtrackVrNvSapFactors.png', 'V_{prnv}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
VrNvSap.XcoordName = Xname
VrNvSap.SimCmap = 'bwr'
VrNvSap.CreateFig()

# Ptrack Vt, vortex removed, press coords
SimCspecs = [ -10, 20, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -4, 4, 11 ]
Xname = '/xl_coords'

VtNvPressPsap = ffc.TrackXsectionPress(VelFname, '/ps_v_nv_lite_p', 'Plots.py/FsFigPtrackVtNvPressPsapFactors.png', 'V_{ptnv}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
VtNvPressPsap.XcoordName = Xname
VtNvPressPsap.CreateFig()

VtNvPressSap = ffc.TrackXsectionPress(VelFname, '/s_v_nv_lite_p', 'Plots.py/FsFigPtrackVtNvPressSapFactors.png', 'V_{ptnv}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
VtNvPressSap.XcoordName = Xname
VtNvPressSap.CreateFig()

# Ptrack Vr, vortex removed
SimCspecs = [ -5, 5, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -2, 2, 11 ]

VrNvPressPsap = ffc.TrackXsectionPress(VelFname, '/ps_u_nv_lite_p', 'Plots.py/FsFigPtrackVrNvPressPsapFactors.png', 'V_{prnv}', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
VrNvPressPsap.XcoordName = Xname
VrNvPressPsap.SimCmap = 'bwr'
VrNvPressPsap.CreateFig()

VrNvPressSap = ffc.TrackXsectionPress(VelFname, '/s_u_nv_lite_p', 'Plots.py/FsFigPtrackVrNvPressSapFactors.png', 'V_{prnv}', 'SAP', SimCspecs, FactCspecs, 'ptrack')
VrNvPressSap.XcoordName = Xname
VrNvPressSap.SimCmap = 'bwr'
VrNvPressSap.CreateFig()





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

## Temp C
#SimCspecs = [ -50, 30, 11 ]
#FactCspecs = [ -3, 3, 11 ]
#Ylim = [ 0, 10 ]
#Yticks = [ 0, 2, 4, 6, 8, 10 ]
#
#TempcPsap = ffc.TrackXsection(AvgFname, '/ps_tempc', 'Plots.py/FsFigPtrackTempcPsapFactors.png', 'T,\\ ^{\\circ}C', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
#TempcPsap.Ylim = Ylim
#TempcPsap.Yticks = Yticks
#TempcPsap.CreateFig()
#
#TempcSap = ffc.TrackXsection(AvgFname, '/s_tempc', 'Plots.py/FsFigPtrackTempcSapFactors.png', 'T,\\ ^{\\circ}C', 'SAP', SimCspecs, FactCspecs, 'ptrack')
#TempcSap.Ylim = Ylim
#TempcSap.Yticks = Yticks
#TempcSap.CreateFig()

## Temp C, pressure coords
#SimCspecs = [ -50, 30, 11 ]
#FactCspecs = [ -3, 3, 11 ]
#
#TempcPressPsap = ffc.TrackXsectionPress(AvgFname, '/ps_tempc_p', 'Plots.py/FsFigPtrackTempcPressPsapFactors.png', 'T,\\ ^{\\circ}C', 'PSAP', SimCspecs, FactCspecs, 'ptrack')
#TempcPressPsap.CreateFig()
#
#TempcPressSap = ffc.TrackXsectionPress(AvgFname, '/s_tempc_p', 'Plots.py/FsFigPtrackTempcPressSapFactors.png', 'T,\\ ^{\\circ}C', 'SAP', SimCspecs, FactCspecs, 'ptrack')
#TempcPressSap.CreateFig()

