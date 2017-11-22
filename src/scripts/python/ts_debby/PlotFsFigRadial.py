#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import FsFigLine as ffl


# Tangential wind speed
Fname = "DIAGS/hist_meas_az_speed_<SIM>.h5"
Vname = "/all_ps_speed_t_maxlev"
Scale = 1.0
Offset = 0.0
Label = r'$V_t\ (ms^{-1})$'
SimLimits = [ 0, 20 ]
FacLimits = [ -2, 2 ]
OutFname = "Plots.py/FsFigRadialVtPsapFactors.png"
Ptitle = r'$V_t$'

Lplot = ffl.Radial(Fname, Vname, Scale, Offset, Label, SimLimits, FacLimits, OutFname, Ptitle)
Lplot.CreateFig()

Fname = "DIAGS/hist_meas_az_speed_<SIM>.h5"
Vname = "/all_s_speed_t_maxlev"
Scale = 1.0
Offset = 0.0
Label = r'$V_t\ (ms^{-1})$'
SimLimits = [ 0, 20 ]
FacLimits = [ -3, 3 ]
OutFname = "Plots.py/FsFigRadialVtSapFactors.png"
Ptitle = r'$V_t$'

Lplot = ffl.Radial(Fname, Vname, Scale, Offset, Label, SimLimits, FacLimits, OutFname, Ptitle)
Lplot.CreateFig()

