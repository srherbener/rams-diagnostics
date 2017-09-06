#!/usr/bin/env python3
import numpy as np

def SetTimeString():
    # Set TimeString to the beginning of the TS Debby spin up simulation
    TimeString = '2012-01-01 00:00:00 00:00'
    return TimeString

def SetColorScheme():
    ColorScheme = {
        'RCE_1km'          : 'black',
        'RCE_1km_SM'       : 'cyan',
        'RCE_1km_DM'       : 'magenta',
        'RCE_1km_DP'       : 'sandybrown',
        'RCE_1km_IPL2'     : 'black',
        'RCE_1km_SM_IPL2'  : 'cyan',
        'RCE_1km_DM_IPL2'  : 'magenta',
        'RCE_1km_DP_IPL2'  : 'sandybrown',
        'RCE_3km_1mom'     : 'black',
        'RCE_3km_2mom'     : 'cyan',
        'RCE_3km_2mom_dm'  : 'magenta',
        'RCE_3km_2mom_db'  : 'sandybrown',
        }
    return ColorScheme

def SetLabelScheme():
    LabelScheme = {
        'RCE_1km'          : '2M',
        'RCE_1km_SM'       : '1M',
        'RCE_1km_DM'       : '2M, DM',
        'RCE_1km_DP'       : '2M, DB',
        'RCE_1km_IPL2'     : '2M',
        'RCE_1km_SM_IPL2'  : '1M',
        'RCE_1km_DM_IPL2'  : '2M, DM',
        'RCE_1km_DP_IPL2'  : '2M, DB',
        'RCE_3km_1mom'     : '1M',
        'RCE_3km_2mom'     : '2M',
        'RCE_3km_2mom_dm'  : '2M, DM',
        'RCE_3km_2mom_db'  : '2M, DB',
        }
    return LabelScheme

def SetHcoords(Min, Max, N):
    Val = Min
    Inc = (Max - Min) / float(N - 1)
    Coords = np.zeros(N)
    for i in np.arange(N):
        Coords[i] = Val
        Val = Val + Inc

    return Coords

def SetZcoords():
    Zlevels = [
        -7.3,    10.4,    47.7,   100.9,   156.5,   218.6,   287.8,   364.9,   450.8,   546.5,
       653.0,   771.7,   903.9,  1051.2,  1215.3,  1398.0,  1601.6,  1828.4,  2081.0,  2362.4,
      2675.8,  3024.9,  3413.9,  3847.6,  4329.1,  4831.9,  5331.9,  5831.9,  6331.9,  6831.9,
      7331.9,  7831.9,  8331.9,  8831.9,  9331.9,  9831.9, 10331.9, 10831.9, 11331.9, 11831.9,
     12331.9, 12831.9, 13331.9, 13831.9, 14331.9, 14831.9, 15331.9, 15831.9, 16331.9, 16831.9,
     17331.9, 17831.9, 18331.9, 18831.9, 19331.9, 19831.4, 20332.4, 20862.9, 21449.2, 22094.1,
     22803.5, 23583.9, 24442.3, 25386.5,
     ]
    return Zlevels
