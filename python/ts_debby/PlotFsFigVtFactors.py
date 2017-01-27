#!/usr/bin/env python3

import os
import sys
sys.path.append("{0:s}/etc/python/common".format(os.environ['HOME']))
sys.path.append("{0:s}/etc/python/ts_debby".format(os.environ['HOME']))

import matplotlib.pyplot as plt
import numpy as np
import h5py
import PlotUtils as plu
import ConfigTsd as conf

LabelScheme = conf.SetLabelScheme()

CaseList = [ 
    'TSD_NONSAL_NODUST',
    'TSD_SAL_NODUST',
    'TSD_NONSAL_DUST',
    'TSD_SAL_DUST'
    ]
Ncases = len(CaseList)

FacList = [
    'FAC_SAL',
    'FAC_DUST',
    'FAC_INT',
    ]
Nfacs = len(FacList)

# Read in Vt data
VtTitles = [ ]
for i in range(0, Ncases):
    Case  = CaseList[i]

    # Get the PSAP data
    InFname = "DIAGS/storm_xsections_{0:s}.h5".format(Case)
    InVname = "/all_ps_speed_t"
    print("Reading {0:s} ({1:s})".format(InFname, InVname))

    InFile = h5py.File(InFname, mode='r')

    # grab the coordinate values, and initialize the wind array
    if (i == 0):
        X = InFile['/x_coords'][:] / 1000.0  # radius in km
        Z = InFile['/z_coords'][:] / 1000.0  # height in km
        Nx = len(X)
        Nz = len(Z)
        VtPsap = np.zeros( [ Nz, Nx, Ncases ], dtype=np.float_)
        VtSap  = np.zeros( [ Nz, Nx, Ncases ], dtype=np.float_)

    VtPsap[ :, :, i ] = InFile[InVname][:,:]
    
    # Get the SAP data
    InVname = "/all_s_speed_t"
    print("Reading {0:s} ({1:s})".format(InFname, InVname))
    VtSap[ :, :, i ] = InFile[InVname][:,:]

    InFile.close()

    VtTitles.append(r'$V_t ({0:s})$'.format(LabelScheme[Case]))

print("")

VtPmarkers = [ 'a', 'c', 'e', 'g' ]

# Build the factors from the input Vt cross sections
#   The third dimension in VtPsap is the case which
#   is mapped as:
#
#      index     case
#        0       NSND
#        1        SND
#        2        NSD
#        3         SD
#
# Factors are in this order:
#
#      index     factor            formula
#        0       FAC_SAL     SND - NSND
#        1      FAC_DUST     NSD - NSND
#        2       FAC_INT     SD - (SND + NSD) + NSND

FacTitles = [ ]
FacPsap = np.zeros( [ Nz, Nx, Nfacs ], dtype=np.float_)
FacSap  = np.zeros( [ Nz, Nx, Nfacs ], dtype=np.float_)
for i in range(0, Nfacs):
    Factor = FacList[i]

    if (Factor == 'FAC_SAL'):
        # SND - NSND
        FacPsap[:,:,i] = VtPsap[:,:,1] - VtPsap[:,:,0]
        FacSap[:,:,i] = VtSap[:,:,1] - VtSap[:,:,0]
    elif (Factor == 'FAC_DUST'):
        # NSD - NSND
        FacPsap[:,:,i] = VtPsap[:,:,2] - VtPsap[:,:,0]
        FacSap[:,:,i] = VtSap[:,:,2] - VtSap[:,:,0]
    elif (Factor == 'FAC_INT'):
        # SD - (SND + NSD) + NSND
        FacPsap[:,:,i] = VtPsap[:,:,3] - (VtPsap[:,:,1] + VtPsap[:,:,2]) + VtPsap[:,:,0]
        FacSap[:,:,i] = VtSap[:,:,3] - (VtSap[:,:,1] + VtSap[:,:,2]) + VtSap[:,:,0]
    else:
        print("ERROR: PlotFsFigVtFactors: Attempted to calculate too many factors")
        print("ERROR: PlotFsFigVtFactors: Stopping!")
        sys.exit(1)

    FacTitles.append(r'$V_t ({0:s})$'.format(LabelScheme[Factor]))

FacPmarkers = [ 'd', 'f', 'h' ]

# 7 panel plot
#   Vt on the left side, 4 panels
#   Factors on right side, 3 panels
FigPsap = plt.figure()
FigSap = plt.figure()

Pdir = 'Plots.py'

# Place axes for all panels
AxW = 0.40
AxH = 0.16
AxLeftLlx = 0.07
AxRightLlx = 0.53
AxLly = [ 0.10, 0.32, 0.54, 0.76 ]

VtAxPsap = [ ]
VtAxPsap.append(FigPsap.add_axes([ AxLeftLlx, AxLly[3], AxW, AxH ]))
VtAxPsap.append(FigPsap.add_axes([ AxLeftLlx, AxLly[2], AxW, AxH ]))
VtAxPsap.append(FigPsap.add_axes([ AxLeftLlx, AxLly[1], AxW, AxH ]))
VtAxPsap.append(FigPsap.add_axes([ AxLeftLlx, AxLly[0], AxW, AxH ]))

VtAxSap = [ ]
VtAxSap.append(FigSap.add_axes([ AxLeftLlx, AxLly[3], AxW, AxH ]))
VtAxSap.append(FigSap.add_axes([ AxLeftLlx, AxLly[2], AxW, AxH ]))
VtAxSap.append(FigSap.add_axes([ AxLeftLlx, AxLly[1], AxW, AxH ]))
VtAxSap.append(FigSap.add_axes([ AxLeftLlx, AxLly[0], AxW, AxH ]))

FacAxPsap = [ ]
FacAxPsap.append(FigPsap.add_axes([ AxRightLlx, AxLly[2], AxW, AxH ]))
FacAxPsap.append(FigPsap.add_axes([ AxRightLlx, AxLly[1], AxW, AxH ]))
FacAxPsap.append(FigPsap.add_axes([ AxRightLlx, AxLly[0], AxW, AxH ]))

FacAxSap = [ ]
FacAxSap.append(FigSap.add_axes([ AxRightLlx, AxLly[2], AxW, AxH ]))
FacAxSap.append(FigSap.add_axes([ AxRightLlx, AxLly[1], AxW, AxH ]))
FacAxSap.append(FigSap.add_axes([ AxRightLlx, AxLly[0], AxW, AxH ]))

Cfilled = 1  # 0 - no fill, 1 - fill

Xaxis = plu.AxisConfig('x', [ 0, 500 ], "Radius (km)")
Xaxis.fontsize = 12

Yaxis = plu.AxisConfig('y', [ 0, 6 ], "Z (km)")
Yaxis.ticks = [ 0, 2, 4, 6 ]
Yaxis.fontsize = 12

# Vt cross sections
Cmap = 'nipy_spectral'
Cmin = 0
Cmax = 20
Cnum = 21
Clim = [ Cmin, Cmax ]
Clevs = np.linspace(Cmin, Cmax, num=Cnum)
Xshow = [ 0, 0, 0, 1 ]
Yshow = [ 1, 1, 1, 1 ]
for i in range(0, Ncases):
    Ptitle = plu.TitleConfig(VtPmarkers[i], VtTitles[i])
    Ptitle.fontsize = 18

    Xaxis.show = Xshow[i]
    Yaxis.show = Yshow[i]
    plu.PlotContour(VtAxPsap[i], X, Z, VtPsap[:,:,i], Ptitle, Xaxis, Yaxis, Cmap, Clim, Clevs, Cfilled)
    plu.PlotContour(VtAxSap[i], X, Z, VtSap[:,:,i], Ptitle, Xaxis, Yaxis, Cmap, Clim, Clevs, Cfilled)

# factor cross sections
Cmap = 'bwr'
Cmin = -4
Cmax = 4
Cnum = 21
Clim = [ Cmin, Cmax ]
Clevs = np.linspace(Cmin, Cmax, num=Cnum)
Xshow = [ 0, 0, 1 ]
Yshow = [ 0, 0, 0 ]
for i in range(0, Nfacs):
    Ptitle = plu.TitleConfig(FacPmarkers[i], FacTitles[i])
    Ptitle.fontsize = 18

    Xaxis.show = Xshow[i]
    Yaxis.show = Yshow[i]
    plu.PlotContour(FacAxPsap[i], X, Z, FacPsap[:,:,i], Ptitle, Xaxis, Yaxis, Cmap, Clim, Clevs, Cfilled)
    plu.PlotContour(FacAxSap[i], X, Z, FacSap[:,:,i], Ptitle, Xaxis, Yaxis, Cmap, Clim, Clevs, Cfilled)



OutFile = "{0:s}/FsFigVtPsapFactors.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
FigPsap.savefig(OutFile)

OutFile = "{0:s}/FsFigVtSapFactors.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
FigSap.savefig(OutFile)

