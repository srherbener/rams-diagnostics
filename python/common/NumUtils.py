#!/usr/bin/env python3
#
# Numeric utilities
#

import numpy as np
from scipy.signal import filtfilt

#########################################################################################
# SmoothLine
#
# Remove noise from a linear data series.
#
def SmoothLine(Var):
    # Use filtfilt which goes forward once and backward once in order to eliminate
    # phase shift. The Gustafsson method is used for extending the endpoints of the
    # series in order to minimize distortion at the ends of the smoothed line.
    Flen = 5

    a = 1.0                           # coeffs for a 'box' filter
    b = np.ones(Flen) / float(Flen)

    VarSmooth = filtfilt(b, a, Var, method="gust")

    return VarSmooth

