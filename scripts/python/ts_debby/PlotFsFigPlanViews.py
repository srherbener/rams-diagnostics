#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import FsFigContour as ffc


InFname = 'DIAGS/wind_shear_nv_lite_<SIM>.h5'

# Storm locations based on TSD_SAL_DUST sim
# Entries are for t = 10,20,30,40,50,60 hours
StormLocs = [
    [ -25.9, 14.5 ],
    [ -27.8, 16.1 ],
    [ -30.3, 16.5 ],
    [ -32.8, 17.3 ],
    [ -35.3, 18.5 ],
    [ -37.8, 19.3 ],
    ]

# Wind shear magnitude
SimCspecs = [ 0, 20, 11 ]  # color specs: [ Cmin, Cmax, Cnum ], Cnum is number of contour levels
FactCspecs = [ -2, 2, 11 ]

MagShear10 = ffc.PlanView(InFname, '/mag_shear_10', 'Plots.py/FsFigWindShear10.png', '|Shear|', '22Aug, 16Z', SimCspecs, FactCspecs, StormLoc=StormLocs[0])
MagShear10.CreateFig()

MagShear20 = ffc.PlanView(InFname, '/mag_shear_20', 'Plots.py/FsFigWindShear20.png', '|Shear|', '23Aug, 02Z', SimCspecs, FactCspecs, StormLoc=StormLocs[1])
MagShear20.CreateFig()

MagShear30 = ffc.PlanView(InFname, '/mag_shear_30', 'Plots.py/FsFigWindShear30.png', '|Shear|', '23Aug, 12Z', SimCspecs, FactCspecs, StormLoc=StormLocs[2])
MagShear30.CreateFig()

MagShear40 = ffc.PlanView(InFname, '/mag_shear_40', 'Plots.py/FsFigWindShear40.png', '|Shear|', '23Aug, 22Z', SimCspecs, FactCspecs, StormLoc=StormLocs[3])
MagShear40.CreateFig()

MagShear50 = ffc.PlanView(InFname, '/mag_shear_50', 'Plots.py/FsFigWindShear50.png', '|Shear|', '24Aug, 08Z', SimCspecs, FactCspecs, StormLoc=StormLocs[4])
MagShear50.CreateFig()

MagShear60 = ffc.PlanView(InFname, '/mag_shear_60', 'Plots.py/FsFigWindShear60.png', '|Shear|', '24Aug, 18Z', SimCspecs, FactCspecs, StormLoc=StormLocs[5])
MagShear60.CreateFig()
