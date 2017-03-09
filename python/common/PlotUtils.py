##############################################################
# Plotting utilities
#
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys

# classes
class AxisConfig:
    '''Class for configuring plot axes'''
    def __init__(self, axid, lim, label):
        self.axid = axid
        self.lim = lim
        self.label = label
        self.scale = "linear"
        self.show = 1
        self.ticks = [ ]
        self.ticklabels = [ ]
        self.fontsize = 8

    ######################################################
    # Method to configure axis according to class variable
    # settings.
    def config(self, Paxes):
        if (self.lim != [ ]):
            if (self.axid == 'x'):
                Paxes.set_xlim(self.lim)
            elif (self.axid == 'y'):
                Paxes.set_ylim(self.lim)

        if (self.label != ""):
            if (self.axid == 'x'):
                Paxes.set_xlabel(self.label, fontsize=self.fontsize)
            elif (self.axid == 'y'):
                Paxes.set_ylabel(self.label, fontsize=self.fontsize)

        if (self.scale != ""):
            if (self.axid == 'x'):
                Paxes.set_xscale(self.scale)
            elif (self.axid == 'y'):
                Paxes.set_yscale(self.scale)

        if (self.ticks != [ ]):
            if (self.axid == 'x'):
                Paxes.set_xticks(self.ticks)
            elif (self.axid == 'y'):
                Paxes.set_yticks(self.ticks)

        if (self.ticklabels != [ ]):
            if (self.axid == 'x'):
                Paxes.set_xticklabels(self.ticklabels)
            elif (self.axid == 'y'):
                Paxes.set_yticklabels(self.ticklabels)

        # Won't necessarily be explicitly setting the xtick or ytick labels
        # so make sure the default labels have their font size set
        if (self.fontsize > 0):
            if (self.axid == 'x'):
                for l in Paxes.get_xticklabels():
                    l.set_fontsize(self.fontsize)
            elif (self.axid == 'y'):
                for l in Paxes.get_yticklabels():
                    l.set_fontsize(self.fontsize)

        if (self.show == 0):
            if (self.axid == 'x'):
                plt.setp(Paxes.axes.get_xticklabels(), visible=False)
                Paxes.set_xlabel('')
            elif (self.axid == 'y'):
                plt.setp(Paxes.axes.get_yticklabels(), visible=False)
                Paxes.set_ylabel('')


class TitleConfig:
    '''Class for configuring plot title'''
    def __init__(self, pmarker, title):
       self.pmarker = pmarker
       self.title = title
       self.fontsize = 14

    #######################################################
    # Method to set title. If a panel marker does not exist,
    # the title will be centered, otherwise it will be
    # left justified.
    def set(self, Paxes):
        if (self.pmarker == ''):
            Paxes.set_title(self.title, fontsize=self.fontsize)
        else:
            Paxes.set_title("({0:s}) {1:s}".format(self.pmarker, self.title), loc='left', fontsize=self.fontsize)


class LegendConfig:
    '''Class for configuring plot legend'''
    def __init__(self, text, loc):
        self.text = text
        self.loc = loc
        self.fontsize = 8
        self.ncol = 1
        self.bbox = [ ]

    #######################################################
    # Method to set the legend. If the location is 'none',
    # then no legend is set.
    def set(self, Paxes):
        if (self.loc != 'none'):
            if (self.bbox == [ ]):
                Paxes.legend(loc=self.loc, ncol=self.ncol)
            else:
                Paxes.legend(bbox_to_anchor=self.bbox, loc=self.loc, ncol=self.ncol)
            
            if (self.fontsize > 0):
                for i in Paxes.get_legend().get_texts():
                    i.set_fontsize(self.fontsize)

class ContourConfig:
    '''Class for configuring contouring'''
    def __init__(self, cmin, cmax, cnum, cmap, cfilled, ctype):
        self.cmin    = cmin
        self.cmax    = cmax
        self.cnum    = cnum
        self.cmap    = cmap
        self.cfilled = cfilled
        self.ctype   = ctype

        if (self.ctype == 'linear'):
            self.clevels = np.linspace(self.cmin, self.cmax, num=self.cnum)
        elif (self.ctype == 'log'):
            self.clevels = np.logspace(self.cmin, self.cmax, num=self.cnum)
        else:
            print("Error: ContourConfig: contour type (ctype) must be one of 'log', 'linear': {0:s}".format(self.ctype))
            sys.exit(1)


##############################################################
# PlotLine()
#
# This routine will plot a set of lines on the same axes.
#
def PlotLine(Paxes, X, Y, Ptitle, Xaxis, Yaxis, Legend, Colors):

    # X and Y need to contain line data in each column. Therefore, the length of the columns
    # (ie, the number of rows) of X and Y need to match. X needs to be a vector, Y can either
    # be a vector or an array.

    LineW = 2

    if (X.ndim == 1):
        Xnpts = X.size
        Xnlines = 1
    else:
        print("ERROR: PlotLine: X needs to be a vector")
        print("ERROR: PlotLine: skipping this plot")
        print("")
        return
    if (Y.ndim == 1):
        Ynpts = Y.size
        Ynlines = 1
    else:
        ( Ynpts, Ynlines ) = Y.shape

    if (Xnpts != Ynpts):
        print("ERROR: PlotLine: number of rows of X and Y need to match")
        print("ERROR: PlotLine: skipping this plot")
        print("")
        return

    # Start drawing
    if (Ynlines == 1):
        Paxes.plot(X, Y, color=Colors[0], linewidth=LineW, label=Legend.text[0])
    else:
        for i in range(0, Ynlines):
            Paxes.plot(X, Y[:,i], color=Colors[i], linewidth=LineW, label=Legend.text[i])

    Ptitle.set(Paxes)
    Xaxis.config(Paxes)
    Yaxis.config(Paxes)
    Legend.set(Paxes)


#####################################################################
# PlotSplitBgraph()
#
# This routine will plot a bar graph where each entry is a separate
# bar on the graph. One use for this would be to be able to depict
# accumulation of factors where individual factors can be positive
# or negative (doing this as stacked bars will cause factors to
# overlap messing up the stack for example).
#
def PlotSplitBgraph(Paxes, Xbars, Ybars1, Ybars2, Ptitle, Xaxis, Yaxis, Legend, Colors):

    # Y1 and Y2 contain the start and end heights of the bars, respectively.
    Nbars = Xbars.size
    DeltaX = (Xbars[1] - Xbars[0]) * 0.4  # assume equal spacing for now

    for i in range(0, Nbars):
        X1 = Xbars[i] - DeltaX
        X2 = Xbars[i] + DeltaX
        Y1 = Ybars1[i]
        Y2 = Ybars2[i]

        Width = X2 - X1
        Height = Y2 - Y1
        Paxes.add_patch(pat.Rectangle([ X1, Y1 ], Width, Height, color=Colors[i]))
 
    Ptitle.set(Paxes)
    Xaxis.config(Paxes)
    Yaxis.config(Paxes)
    Legend.set(Paxes)

#####################################################################
# PlotContour()
#
# This routine will plot a 2D contour plot.
#
def PlotContour(Paxes, X, Y, Z, Ptitle, Xaxis, Yaxis, Cspecs):

    # Make the contour plot
    if (Cspecs.cfilled):
        cplot = Paxes.contourf(X, Y, Z, linestyle='none', extend='both', levels=Cspecs.clevels,
            vmin=Cspecs.cmin, vmax=Cspecs.cmax, cmap=Cspecs.cmap)
    else:
        cplot = Paxes.contour(X, Y, Z, extend='both', levels=Cspecs.clevels,
            vmin=Cspecs.cmin, vmax=Cspecs.cmax, cmap=Cspecs.cmap)

    plt.colorbar(cplot, ax=Paxes, aspect=10)

    Ptitle.set(Paxes)
    Xaxis.config(Paxes)
    Yaxis.config(Paxes)

