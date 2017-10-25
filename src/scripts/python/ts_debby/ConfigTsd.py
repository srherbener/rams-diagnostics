#!/usr/bin/env python3

def SetTimeString():
    # Set TimeString to the beginning of the TS Debby spin up simulation
    TimeString = '2006-08-20 12:00:00 00:00'
    return TimeString

def SetRhoAir():
    RhoAir = [
        1.196, 1.191, 1.185, 1.179, 1.172, 1.165, 1.157, 1.147, 1.136, 1.124,
        1.111, 1.097, 1.082, 1.066, 1.051, 1.034, 1.017, 1.000, 0.982, 0.963,
        0.945, 0.925, 0.905, 0.883, 0.860, 0.833, 0.804, 0.773, 0.741, 0.711,
        0.679, 0.647, 0.612, 0.577, 0.541, 0.505, 0.469, 0.432, 0.394, 0.355,
        0.316, 0.279, 0.244, 0.210, 0.179, 0.150, 0.126, 0.105, 0.087, 0.073,
        0.062, 0.052, 0.044, 0.038, 0.032, 0.027
        ]
    return RhoAir

def SetColorScheme():
    ColorScheme = {
        'TSD_SAL_DUST'     : 'black',
        'TSD_SAL_NODUST'   : 'cyan',
        'TSD_NONSAL_DUST'  : 'magenta',
        'TSD_NONSAL_NODUST': 'sandybrown',
        'FAC_SAL'          : 'cyan',
        'FAC_DUST'         : 'magenta',
        'FAC_INT'          : 'seagreen',
        }
    return ColorScheme

def SetLabelScheme():
    LabelScheme = {
        'TSD_SAL_DUST'     : 'SD',
        'TSD_SAL_NODUST'   : 'SND',
        'TSD_NONSAL_DUST'  : 'NSD',
        'TSD_NONSAL_NODUST': 'NSND',
        'FAC_SAL'          : 'FS',
        'FAC_DUST'         : 'FD',
        'FAC_INT'          : 'FSD',
        }
    return LabelScheme
