#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import FsFigLine as ffl


# Precip Rate
Fname = "DIAGS/hda_meas_ts_pcprate_<SIM>.h5"
Vname = "/nw_pcprate_ts"
Scale = 1.0
Offset = 0.0
Label = r'$PR\ (mm\ hr^{-1})$'
SimLimits = [ 0, 0.5 ]
FacLimits = [ -0.3, 0.3 ]
OutFname = "Plots.py/FsFigTseriesNwPrecipRateFactors.png"
Ptitle = r'$Precip\ Rate,\ NW$'

Lplot = ffl.TimeSeries(Fname, Vname, Scale, Offset, Label, SimLimits, FacLimits, OutFname, Ptitle)
Lplot.CreateFig()


Fname = "DIAGS/hda_meas_ts_pcprate_<SIM>.h5"
Vname = "/nw_nstorm_pcprate_ts"
Scale = 1.0
Offset = 0.0
Label = r'$PR\ (mm\ hr^{-1})$'
SimLimits = [ 0, 0.3 ]
FacLimits = [ -0.2, 0.2 ]
OutFname = "Plots.py/FsFigTseriesNwNstormPrecipRateFactors.png"
Ptitle = r'$Precip\ Rate,\ NW\ NSTORM$'

Lplot = ffl.TimeSeries(Fname, Vname, Scale, Offset, Label, SimLimits, FacLimits, OutFname, Ptitle)
Lplot.CreateFig()


# Total precip
Fname = "DIAGS/hda_meas_ts_precip_<SIM>.h5"
Vname = "/nw_sum_accpcp_mass"
Scale = 1e-4
Offset = 0.0
Label = r'$Precip\ (m\ X\ 10)$'
SimLimits = [ 0, 200 ]
FacLimits = [ -50, 50 ]
OutFname = "Plots.py/FsFigTseriesNwAccumPrecipFactors.png"
Ptitle = r'$Accum\ Precip,\ NW$'

Lplot = ffl.TimeSeries(Fname, Vname, Scale, Offset, Label, SimLimits, FacLimits, OutFname, Ptitle)
Lplot.CreateFig()


Fname = "DIAGS/hda_meas_ts_precip_<SIM>.h5"
Vname = "/nw_nstorm_sum_accpcp_mass"
Scale = 1e-4
Offset = 0.0
Label = r'$Precip\ (m\ X\ 10)$'
SimLimits = [ 0, 200 ]
FacLimits = [ -50, 50 ]
OutFname = "Plots.py/FsFigTseriesNwNstormAccumPrecipFactors.png"
Ptitle = r'$Accum\ Precip,\ NW\ NSTORM$'

Lplot = ffl.TimeSeries(Fname, Vname, Scale, Offset, Label, SimLimits, FacLimits, OutFname, Ptitle)
Lplot.CreateFig()


# Precipitable water (vertically integrated vapor)
Fname = "DIAGS/hda_meas_ts_vint_vapor_<SIM>.h5"
Vname = "/nw_sum_vint_vapor"
Scale = 1e-6
Offset = 0.0
Label = r'$PW\ (km)$'
SimLimits = [ 3, 8 ]
FacLimits = [ -1, 1 ]
OutFname = "Plots.py/FsFigTseriesNwVintVaporFactors.png"
Ptitle = r'$Precipitable\ Water,\ NW$'

Lplot = ffl.TimeSeries(Fname, Vname, Scale, Offset, Label, SimLimits, FacLimits, OutFname, Ptitle)
Lplot.CreateFig()


Fname = "DIAGS/hda_meas_ts_vint_vapor_<SIM>.h5"
Vname = "/nw_nstorm_sum_vint_vapor"
Scale = 1e-6
Offset = 0.0
Label = r'$Precip\ (km)$'
SimLimits = [ 3, 8 ]
FacLimits = [ -1, 1 ]
OutFname = "Plots.py/FsFigTseriesNwNstormVintVaporFactors.png"
Ptitle = r'$Preciptable\ Water,\ NW\ NSTORM$'

Lplot = ffl.TimeSeries(Fname, Vname, Scale, Offset, Label, SimLimits, FacLimits, OutFname, Ptitle)
Lplot.CreateFig()

