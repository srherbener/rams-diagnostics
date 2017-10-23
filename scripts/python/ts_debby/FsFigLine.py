###################################################################
# Templates for creating line plots that show simulation
# results on the top panel, and factors on the bottom panel.
#

import matplotlib.pyplot as plt
import numpy as np
import h5py
import PlotUtils as plu
import ConfigTsd as conf

# base class - define the panel layout and create the plot
class FsLine:
    '''Class for plotting lines showing sims on the top and factors on the bottom'''

    def __init__(self, InFname, Xname, Yname, OutFname, Ptitle):
        self.InFname = InFname
        self.Xname = Xname
        self.Xscale = 1.0
        self.Xoffset = 0.0
        self.Yname = Yname
        self.Yscale = 1.0
        self.Yoffset = 0.0

        self.OutFname = OutFname

        self.Title = Ptitle

        self.Xlabel = ''
        self.Xmin = 0.0
        self.Xmax = 10.0
        self.Xticks = []
        self.XtickLabels = []
        self.SimXshow = [ 0 ]
        self.FacXshow = [ 1 ]

        self.Ylabel = ''
        self.Ymin = 0.0
        self.Ymax = 10.0
        self.Yticks = []
        self.YtickLabels = []
        self.SimYshow = [ 1 ]
        self.FacYshow = [ 1 ]

        self.AxFsize = 12

        self.SimList = [ 
            'TSD_NONSAL_NODUST',
            'TSD_SAL_NODUST',
            'TSD_NONSAL_DUST',
            'TSD_SAL_DUST'
            ]
        self.Nsims = len(self.SimList)

        self.FacList = [
            'FAC_SAL',
            'FAC_DUST',
            'FAC_INT',
            ]
        self.Nfacs = len(self.FacList)

        self.TitleFsize = 18
        self.SimPmarkers = [ 'a' ]
        self.FacPmarkers = [ 'b' ]

        self.Fig = 0
        self.SimAx = []
        self.FacAx = []

    def CreateFig(self):

        # get labels and colors for sims and factors
        LabelScheme = conf.SetLabelScheme()
        ColorScheme = conf.SetColorScheme()

        # Collect input data in an array: (y,sim)
        # Create factors on the y values, assume all sims have identical x values
        # and all sims have identical length y values
        SimLabels = []
        SimColors = []
        for i in range(self.Nsims):
            Sim = self.SimList[i]

            SimLabels.append(LabelScheme[Sim])
            SimColors.append(ColorScheme[Sim])
        
            InFname = self.InFname.replace("<SIM>", Sim)
            Xname = self.Xname
            Yname = self.Yname
            print("Reading {0:s} ({1:s})".format(InFname, Xname))
            print("Reading {0:s} ({1:s})".format(InFname, Yname))
        
            InFile = h5py.File(InFname, mode='r')
            X = InFile[Xname][...] * self.Xscale + self.Xoffset
            Y = InFile[Yname][...] * self.Yscale + self.Yoffset

            if (i == 0):
                Nx = len(X)
                Ny = len(Y)
                Yvals = np.zeros([ Ny, self.Nsims ], dtype=np.float_)

            Yvals[:,i] = Y
            
            InFile.close()
        
        print("")

        # Build the factors from the input cross sections
        #   The second dimension in Yvals is the sim which
        #   is mapped as:
        #
        #      index     sim
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
        
        FacLabels = []
        FacColors = []
        FacVals = np.zeros( [ Ny, self.Nfacs ], dtype=np.float_)
        for i in range(self.Nfacs):
            Factor = self.FacList[i]

            FacLabels.append(LabelScheme[Factor])
            FacColors.append(ColorScheme[Factor])
        
            if (Factor == 'FAC_SAL'):
                # SND - NSND
                FacVals[:,i] = Yvals[:,1] - Yvals[:,0]
            elif (Factor == 'FAC_DUST'):
                # NSD - NSND
                FacVals[:,i] = Yvals[:,2] - Yvals[:,0]
            elif (Factor == 'FAC_INT'):
                # SD - (SND + NSD) + NSND
                FacVals[:,i] = Yvals[:,3] - (Yvals[:,1] + Yvals[:,2]) + Yvals[:,0]
            else:
                print("ERROR: PlotXsection: Attempted to calculate too many factors")
                print("ERROR: PlotXsection: Stopping!")
                sys.exit(1)
        
        # 2 panel plot
        #   Sims on top
        #   Factors on bottom
        self.Fig = plt.figure()

        # legeng specs
        LegBbox = [ 1.02, 1.0 ]
        LegFsize = 12
        LegNcol = 1
        
        # Place axes for all panels
        AxW = 0.70
        AxH = 0.35
        AxLlx = 0.10
        AxLly = [ 0.10, 0.55 ]
        
        self.SimAx.append(self.Fig.add_axes([ AxLlx, AxLly[1], AxW, AxH ]))
        self.FacAx.append(self.Fig.add_axes([ AxLlx, AxLly[0], AxW, AxH ]))
        
        # create axes
        Xaxis = plu.AxisConfig('x', [ self.Xmin, self.Xmax ], self.Xlabel)
        Xaxis.fontsize = self.AxFsize
        if self.Xticks:
             Xaxis.ticks = self.Xticks
        if self.XtickLabels:
             Xaxis.ticklabels = self.XtickLabels
        
        Yaxis = plu.AxisConfig('y', [ self.Ymin, self.Ymax ], self.Ylabel)
        Yaxis.fontsize = self.AxFsize
        if self.Yticks:
             Yaxis.ticks = self.Yticks
        if self.YtickLabels:
             Yaxis.ticklabels = self.YtickLabels

        # Sim lines
        Ptitle = plu.TitleConfig(self.SimPmarkers[0], "SIMS: {0:s}".format(self.Title))
        Ptitle.fontsize = self.TitleFsize

        Legend = plu.LegendConfig(SimLabels, 'upper left')
        Legend.bbox = LegBbox
        Legend.fontsize = LegFsize
        Legend.ncol = LegNcol
        
        Xaxis.show = self.SimXshow[0]
        Yaxis.show = self.SimYshow[0]
        Yaxis.lim = self.SimLimits
        plu.PlotLine(self.SimAx[0], X, Yvals, Ptitle, Xaxis, Yaxis, Legend, SimColors)
        
        # Factor lines
        Ptitle = plu.TitleConfig(self.FacPmarkers[0], "FACTORS: {0:s}".format(self.Title))
        Ptitle.fontsize = self.TitleFsize

        Legend = plu.LegendConfig(FacLabels, 'upper left')
        Legend.bbox = LegBbox
        Legend.fontsize = LegFsize
        Legend.ncol = LegNcol
        
        Xaxis.show = self.FacXshow[0]
        Yaxis.show = self.FacYshow[0]
        Yaxis.lim = self.FacLimits
        plu.PlotLine(self.FacAx[0], X, FacVals, Ptitle, Xaxis, Yaxis, Legend, FacColors)
        
        print("Writing: {0:s}".format(self.OutFname))
        self.Fig.savefig(self.OutFname)
        plt.close()
        print("")


# class for doing time series
class TimeSeries(FsLine):
    '''Class to create time series, sim/factor figure'''

    def __init__(self, InFname, InVname, InScale, InOffset, InLabel, SimLimits, FacLimits, OutFname, Ptitle):
        FsLine.__init__(self, InFname, "/t_coords", InVname, OutFname, Ptitle)

        # x-axis is time in sim time hours
        self.Xscale = 1.0 / 3600.0
        self.Xoffset = -42.0
        self.Xlabel = ''
        self.Xmin = 0.0
        self.Xmax = 60.0
        self.Xticks = [ 6.0, 18.0, 30.0, 42.0, 54.0 ]
        self.XtickLabels = [ "12Z\n22Aug", "0Z\n23Aug", "12Z\n23Aug", "0Z\n24Aug", "12Z\n24Aug" ]

        # y-axis is height in kilometers
        self.Yscale = InScale
        self.Yoffset = InOffset
        self.Ylabel = InLabel

        self.SimLimits = SimLimits
        self.FacLimits = FacLimits

