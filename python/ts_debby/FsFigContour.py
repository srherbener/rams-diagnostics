###################################################################
# Templates for creating 2D contour plots that show simulation
# results in the left column, and factors in the right column.
#
# The panel layout is 4 rows and 2 columns with seven panels being
# filled in as follows.
#
#     NSND    <blank>
#      SND       FS
#      NSD       FD
#       SD       FSD
#
# All of these share the same axes, so want to show the axes in the
# panels along the left and bottom of the figure. Each entry below
# is <ShowX><ShowY> where '1' means to show the axis, and '0' means
# to not show the axis.
#
#       01
#       01       00
#       01       00
#       11       10
#

import matplotlib.pyplot as plt
import numpy as np
import h5py
import PlotUtils as plu
import ConfigTsd as conf

# base class - define the panel layout and create the plot
class FsContour:
    '''Class for plotting 2D contours showing sims on the right and factors on the left'''

    def __init__(self, InFname, InVname, OutFname, Ptitle, Tag):
        self.InFname = InFname
        self.InVname = InVname
        self.InScale = 1.0
        self.InOffset = 0.0

        self.OutFname = OutFname

        self.Title = Ptitle
        self.Tag = Tag

        self.SimCmin = 0.0
        self.SimCmax = 10.0
        self.SimCnum = 11
        self.SimCmap = 'nipy_spectral'
        self.SimCtype = 'linear'

        self.FacCmin = -2
        self.FacCmax = 2
        self.FacCnum = 11 
        self.FacCmap = 'bwr'
        self.FacCtype = 'linear'

        self.Cfilled = True

        self.XcoordName = '/x_coords'
        self.XcoordScale = 1.0
        self.XcoordOffset = 0.0
        self.Xlabel = ''
        self.Xmin = 0.0
        self.Xmax = 10.0
        self.Xticks = []
        self.XtickLabels = []
        self.SimXshow = [ 0, 0, 0, 1 ]
        self.FacXshow = [ 0, 0, 1 ]

        self.YcoordName = '/y_coords'
        self.YcoordScale = 1.0
        self.YcoordOffset = 0.0
        self.Ylabel = ''
        self.Ymin = 0.0
        self.Ymax = 10.0
        self.Yticks = []
        self.YtickLabels = []
        self.SimYshow = [ 1, 1, 1, 1 ]
        self.FacYshow = [ 0, 0, 0 ]

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
        self.SimPmarkers = [ 'a', 'c', 'e', 'g' ]
        self.FacPmarkers = [ 'd', 'f', 'h' ]

        self.Fig = 0
        self.SimAx = []
        self.FacAx = []

    def CreateFig(self):

        # get labels for sims and factors
        LabelScheme = conf.SetLabelScheme()


        # Collect input data in a 3D array: (x,z,sim) 
        SimTitles = []
        for i in range(self.Nsims):
            Sim = self.SimList[i]
        
            InFname = self.InFname.replace("<SIM>", Sim)
            InVname = self.InVname
            print("Reading {0:s} ({1:s})".format(InFname, InVname))
        
            InFile = h5py.File(InFname, mode='r')
        
            # grab the coordinate values, and initialize the input array
            if (i == 0):
                print("Reading {0:s} ({1:s})".format(InFname, self.XcoordName))
                X = InFile[self.XcoordName][:] * self.XcoordScale + self.XcoordOffset
                print("Reading {0:s} ({1:s})".format(InFname, self.YcoordName))
                Y = InFile[self.YcoordName][:] * self.YcoordScale + self.YcoordOffset
                Nx = len(X)
                Ny = len(Y)
                InVar = np.zeros( [ Ny, Nx, self.Nsims ], dtype=np.float_)
        
            InVar[ :, :, i ] = InFile[InVname][:,:] * self.InScale + self.InOffset
            
            InFile.close()
        
            SimTitles.append(r'${0:s}\ ({1:s})$'.format(self.Title, LabelScheme[Sim]))
        
        # Build the factors from the input cross sections
        #   The third dimension in InVar is the sim which
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
        
        FacTitles = [ ]
        FacVar = np.zeros( [ Ny, Nx, self.Nfacs ], dtype=np.float_)
        for i in range(self.Nfacs):
            Factor = self.FacList[i]
        
            if (Factor == 'FAC_SAL'):
                # SND - NSND
                FacVar[:,:,i] = InVar[:,:,1] - InVar[:,:,0]
            elif (Factor == 'FAC_DUST'):
                # NSD - NSND
                FacVar[:,:,i] = InVar[:,:,2] - InVar[:,:,0]
            elif (Factor == 'FAC_INT'):
                # SD - (SND + NSD) + NSND
                FacVar[:,:,i] = InVar[:,:,3] - (InVar[:,:,1] + InVar[:,:,2]) + InVar[:,:,0]
            else:
                print("ERROR: PlotXsection: Attempted to calculate too many factors")
                print("ERROR: PlotXsection: Stopping!")
                sys.exit(1)
        
            FacTitles.append(r'${0:s}\ ({1:s})$'.format(self.Title, LabelScheme[Factor]))
        
        # 7 panel plot
        #   Sims on the left side, 4 panels
        #   Factors on right side, 3 panels
        self.Fig = plt.figure()
        
        # Place axes for all panels
        AxW = 0.40
        AxH = 0.16
        AxLeftLlx = 0.07
        AxRightLlx = 0.53
        AxLly = [ 0.10, 0.32, 0.54, 0.76 ]
        
        self.SimAx.append(self.Fig.add_axes([ AxLeftLlx, AxLly[3], AxW, AxH ]))
        self.SimAx.append(self.Fig.add_axes([ AxLeftLlx, AxLly[2], AxW, AxH ]))
        self.SimAx.append(self.Fig.add_axes([ AxLeftLlx, AxLly[1], AxW, AxH ]))
        self.SimAx.append(self.Fig.add_axes([ AxLeftLlx, AxLly[0], AxW, AxH ]))
        
        self.FacAx.append(self.Fig.add_axes([ AxRightLlx, AxLly[2], AxW, AxH ]))
        self.FacAx.append(self.Fig.add_axes([ AxRightLlx, AxLly[1], AxW, AxH ]))
        self.FacAx.append(self.Fig.add_axes([ AxRightLlx, AxLly[0], AxW, AxH ]))
        
        # for placing text in upper right (blank) space
        # create axes, turn off the axes, then place the text in that axes
        TxtAx = self.Fig.add_axes([ AxRightLlx, AxLly[3], AxW, AxH ])
        TxtAx.set_axis_off()
        
        # create axes and contour specs objects
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

        SimCspecs = plu.ContourConfig(self.SimCmin, self.SimCmax, self.SimCnum, self.SimCmap, self.Cfilled, self.SimCtype)
        FacCspecs = plu.ContourConfig(self.FacCmin, self.FacCmax, self.FacCnum, self.FacCmap, self.Cfilled, self.FacCtype)

        # Sim cross sections
        for i in range(self.Nsims):
            Ptitle = plu.TitleConfig(self.SimPmarkers[i], SimTitles[i])
            Ptitle.fontsize = self.TitleFsize
        
            Xaxis.show = self.SimXshow[i]
            Yaxis.show = self.SimYshow[i]
            plu.PlotContour(self.SimAx[i], X, Y, InVar[:,:,i], Ptitle, Xaxis, Yaxis, SimCspecs)
        
        # factor cross sections
        for i in range(self.Nfacs):
            Ptitle = plu.TitleConfig(self.FacPmarkers[i], FacTitles[i])
            Ptitle.fontsize = self.TitleFsize
        
            Xaxis.show = self.FacXshow[i]
            Yaxis.show = self.FacYshow[i]
            plu.PlotContour(self.FacAx[i], X, Y, FacVar[:,:,i], Ptitle, Xaxis, Yaxis, FacCspecs)
        
        # Place a tag (text) in the empty panel if requested.
        # transAxes says to use the Axes coordinates which always run from 0 to 1
        if (self.Tag):
            TxtAx.text(1, 1, self.Tag, transform=TxtAx.transAxes,
                fontsize=20, horizontalalignment='right', verticalalignment='bottom')

        # If doing a track (PTRACK or STRACK), place labels on x-axis
        if (hasattr(self, "tracktype")):
            if (self.tracktype == 'ptrack'):
                self.SimAx[3].text(0, -0.1, 'C', transform=self.SimAx[3].transAxes, fontsize=14, color='blue',
                    horizontalalignment='left', verticalalignment='top')
                self.SimAx[3].text(1, -0.1, 'D', transform=self.SimAx[3].transAxes, fontsize=14, color='blue',
                    horizontalalignment='right', verticalalignment='top')
                self.FacAx[2].text(0, -0.1, 'C', transform=self.FacAx[2].transAxes, fontsize=14, color='blue',
                    horizontalalignment='left', verticalalignment='top')
                self.FacAx[2].text(1, -0.1, 'D', transform=self.FacAx[2].transAxes, fontsize=14, color='blue',
                    horizontalalignment='right', verticalalignment='top')
            elif (self.tracktype == 'strack'):
                self.SimAx[3].text(0, -0.1, 'A', transform=self.SimAx[3].transAxes, fontsize=14, color='red',
                    horizontalalignment='left', verticalalignment='top')
                self.SimAx[3].text(1, -0.1, 'B', transform=self.SimAx[3].transAxes, fontsize=14, color='red',
                    horizontalalignment='right', verticalalignment='top')
                self.FacAx[2].text(0, -0.1, 'A', transform=self.FacAx[2].transAxes, fontsize=14, color='red',
                    horizontalalignment='left', verticalalignment='top')
                self.FacAx[2].text(1, -0.1, 'B', transform=self.FacAx[2].transAxes, fontsize=14, color='red',
                    horizontalalignment='right', verticalalignment='top')
        
        # If doing a plan view with storm location, mark the storm location with
        # a large white 'X'.
        if (hasattr(self, "StormX") and hasattr(self, "StormY")):
            for i in range(4):
                self.SimAx[i].plot(self.StormX, self.StormY, 'wo')
            for i in range(3):
                self.FacAx[i].plot(self.StormX, self.StormY, 'co')

        print("Writing: {0:s}".format(self.OutFname))
        self.Fig.savefig(self.OutFname)
        plt.close()
        print("")


# class for doing storm cross sections - radius vs height
class StormXsection(FsContour):
    '''Class to create storm cross section, sim/factor figure'''

    def __init__(self, InFname, InVname, OutFname, Ptitle, Tag, SimCspecs, FacCspecs):
        FsContour.__init__(self, InFname, InVname, OutFname, Ptitle, Tag)

        # x-axis is radius in kilometers
        self.XcoordName = '/x_coords'
        self.XcoordScale = 1.0e-3
        self.XcoordOffset = 0.0
        self.Xlabel = 'Radius (km)'
        self.Xmin = 0.0
        self.Xmax = 500.0
        self.Xticks = [ 100.0, 200.0, 300.0, 400.0, 500.0 ]

        # y-axis is height in kilometers
        self.YcoordName = '/z_coords'
        self.YcoordScale = 1.0e-3
        self.YcoordOffset = 0.0
        self.Ylabel = 'Z (km)'
        self.Ymin = 0.0
        self.Ymax = 8.0
        self.Yticks = [ 0.0, 2.0, 4.0, 6.0, 8.0 ]

        # Contour specs
        self.SimCmin = SimCspecs[0]
        self.SimCmax = SimCspecs[1]
        self.SimCnum = SimCspecs[2]

        self.FacCmin = FacCspecs[0]
        self.FacCmax = FacCspecs[1]
        self.FacCnum = FacCspecs[2]

# class for doing track cross sections - length vs height
class TrackXsection(FsContour):
    '''Class to create track cross section, sim/factor figure'''

    def __init__(self, InFname, InVname, OutFname, Ptitle, Tag, SimCspecs, FacCspecs, TrackType):
        FsContour.__init__(self, InFname, InVname, OutFname, Ptitle, Tag)

        self.tracktype = TrackType

        # x-axis is linear distance in kilometers
        self.XcoordName = '/x_coords'
        self.XcoordScale = 1.0
        self.XcoordOffset = 0.0
        self.Xlabel = 'Linear Distance (km)'
        self.Xmin = 0.0
        if (self.tracktype == 'strack'):
            self.Xmax = 2100.0
        elif (self.tracktype == 'ptrack'):
            self.Xmax = 1800.0
        else:
            self.Xmax = 2000.0
        self.Xticks = [ 500.0, 1000.0, 1500.0 ]

        # y-axis is height in kilometers
        self.YcoordName = '/z_coords'
        self.YcoordScale = 1.0e-3
        self.YcoordOffset = 0.0
        self.Ylabel = 'Z (km)'
        self.Ymin = 0.0
        self.Ymax = 10.0
        self.Yticks = [ 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 ]

        # Contour specs
        self.SimCmin = SimCspecs[0]
        self.SimCmax = SimCspecs[1]
        self.SimCnum = SimCspecs[2]

        self.FacCmin = FacCspecs[0]
        self.FacCmax = FacCspecs[1]
        self.FacCnum = FacCspecs[2]

# class for doing track cross sections - length vs pressure
class TrackXsectionPress(FsContour):
    '''Class to create track cross section, sim/factor figure'''

    def __init__(self, InFname, InVname, OutFname, Ptitle, Tag, SimCspecs, FacCspecs, TrackType):
        FsContour.__init__(self, InFname, InVname, OutFname, Ptitle, Tag)

        self.tracktype = TrackType

        # x-axis is linear distance in kilometers
        self.XcoordName = '/x_coords'
        self.XcoordScale = 1.0
        self.XcoordOffset = 0.0
        self.Xlabel = 'Linear Distance (km)'
        self.Xmin = 0.0
        if (self.tracktype == 'strack'):
            self.Xmax = 2100.0
        elif (self.tracktype == 'ptrack'):
            self.Xmax = 1800.0
        else:
            self.Xmax = 2000.0
        self.Xticks = [ 500.0, 1000.0, 1500.0 ]

        # y-axis is pressure in millibars
        self.YcoordName = '/p_coords'
        self.YcoordScale = 1.0
        self.YcoordOffset = 0.0
        self.Ylabel = 'P (mb)'
        self.Ymin = 1000.0
        self.Ymax = 200.0
        self.Yticks = [ 1000.0, 800.0, 600.0, 400.0, 200.0 ] 

        # Contour specs
        self.SimCmin = SimCspecs[0]
        self.SimCmax = SimCspecs[1]
        self.SimCnum = SimCspecs[2]

        self.FacCmin = FacCspecs[0]
        self.FacCmax = FacCspecs[1]
        self.FacCnum = FacCspecs[2]

# class for doing storm hovmollers - time vs height
class StormHovmoller(FsContour):
    '''Class to create storm hovmoller, sim/factor figure'''

    def __init__(self, InFname, InVname, OutFname, Ptitle, Tag, SimCspecs, FacCspecs):
        FsContour.__init__(self, InFname, InVname, OutFname, Ptitle, Tag)

        self.AxFsize = 10

        # x-axis is time in hours of simulation, start time is 06Z, 22Aug
        self.XcoordName = '/t_coords'
        self.XcoordScale = 1.0 / 3600.0
        self.XcoordOffset = -42.0
        self.Xlabel = ''
        self.Xmin = 0.0
        self.Xmax = 60.0
        self.Xticks = [ 6.0, 18.0, 30.0, 42.0, 54.0 ]
        self.XtickLabels = [ ' 12Z\n22Aug', ' 0Z\n23Aug', ' 12Z\n23Aug', ' 0Z\n24Aug', ' 12Z\n24Aug' ]

        # y-axis is height in kilometers
        self.YcoordName = '/z_coords'
        self.YcoordScale = 1.0e-3
        self.YcoordOffset = 0.0
        self.Ylabel = 'Z (km)'
        self.Ymin = 0.0
        self.Ymax = 8.0
        self.Yticks = [ 0.0, 2.0, 4.0, 6.0, 8.0 ]

        # Contour specs
        self.SimCmin = SimCspecs[0]
        self.SimCmax = SimCspecs[1]
        self.SimCnum = SimCspecs[2]

        self.FacCmin = FacCspecs[0]
        self.FacCmax = FacCspecs[1]
        self.FacCnum = FacCspecs[2]

# class for doing plan views (lon, lat)
class PlanView(FsContour):
    '''Class to create plan view sim/factor figure'''

    def __init__(self, InFname, InVname, OutFname, Ptitle, Tag, SimCspecs, FacCspecs, StormLoc=[]):
        FsContour.__init__(self, InFname, InVname, OutFname, Ptitle, Tag)

        if (StormLoc != []):
            self.StormX = StormLoc[0]
            self.StormY = StormLoc[1]

        # x-axis is longitude in degrees
        self.XcoordName = '/x_coords'
        self.XcoordScale = 1.0
        self.XcoordOffset = 0.0
        self.Xlabel = 'Longitude'
        self.Xmin = -39.5
        self.Xmax = -13.5
        self.Xticks = [ -35, -25, -15 ]
        self.XtickLabels = [ '35W', '25W', '15W' ]

        # y-axis is latitude in degrees
        self.YcoordName = '/y_coords'
        self.YcoordScale = 1.0
        self.YcoordOffset = 0.0
        self.Ylabel = 'Latitude'
        self.Ymin =  7.5
        self.Ymax = 23.5
        self.Yticks = [ 10, 15, 20 ]
        self.YtickLabels = [ '10N', '15N', '20N' ]

        # Contour specs
        self.SimCmin = SimCspecs[0]
        self.SimCmax = SimCspecs[1]
        self.SimCnum = SimCspecs[2]

        self.FacCmin = FacCspecs[0]
        self.FacCmax = FacCspecs[1]
        self.FacCnum = FacCspecs[2]

