#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import PlotTcXsection as ptx


InFname = 'DIAGS/storm_xsections_<CASE>.h5'

# Max Vt
SimCspecs = [ 0, 20, 21 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -4, 4, 21 ]
VtPsap = ptx.StormXsection(InFname, '/all_ps_speed_t', 'V_t', 'PSAP', SimCspecs, FactCspecs, 'Plots.py/FsFigVtPsapFactors.png')
VtPsap.PlotXsection()
VtSap = ptx.StormXsection(InFname, '/all_s_speed_t', 'V_t', 'SAP', SimCspecs, FactCspecs, 'Plots.py/FsFigVtSapFactors.png')
VtSap.PlotXsection()

