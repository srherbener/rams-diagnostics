#!/usr/bin/env python3

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
