##############################################################
# Plotting utilities
#
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.basemap import Basemap

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
    def config_axis(self, Paxes):
        if (self.lim != [ ]):
            if (self.axid == 'x'):
                Paxes.set_xlim(self.lim)
            elif (self.axid == 'y'):
                Paxes.set_ylim(self.lim)

        if (self.label != ""):
            if (self.axid == 'x'):
                Paxes.set_xlabel(self.label)
            elif (self.axid == 'y'):
                Paxes.set_ylabel(self.label)

        if (self.scale != ""):
            if (self.axid == 'x'):
                Paxes.set_xscale(self.scale)
            elif (self.axid == 'y'):
                Paxes.set_yscale(self.scale)


class TitleConfig:
    '''Class for configuring plot title'''
    def __init__(self, pmarker, title):
       self.pmarker = pmarker
       self.title = title

    #######################################################
    # Method to left justify title based on existence of
    # a panel marker.
    def leftjust(self, Paxes):
        if (self.pmarker == ''):
            Paxes.set_title(self.title)
        else:
            Paxes.set_title("({0:s}) {1:s}".format(self.pmarker, self.title), loc='left')


##############################################################
# PlotLine()
#
# This routine will plot a set of lines on the same axes.
#
def PlotLine(Paxes, X, Y, Ptitle, Xaxis, Yaxis, Fsize, LegText, LegLoc, Colors):

    # X and Y need to contain line data in each column. Therefore, the length of the columns
    # (ie, the number of rows) of X and Y need to match. X needs to be a vector, Y can either
    # be a vector or an array.

    TitleFsize = Fsize
    LabelFsize = Fsize * 0.8
    LegendFsize = Fsize * 0.6

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

    # Set the font sizes
    FontParams = { 'legend.fontsize': LegendFsize,
                   'axes.titlesize': TitleFsize,
                   'axes.labelsize': LabelFsize,
                   'xtick.labelsize': LabelFsize,
                   'ytick.labelsize': LabelFsize }
    pylab.rcParams.update(FontParams)

    # Start drawing
    if (Ynlines == 1):
        Paxes.plot(X, Y, color=Colors[0], linewidth=LineW, label=LegText[0])
    else:
        for i in range(0, Ynlines):
            Paxes.plot(X, Y[:,i], color=Colors[i], linewidth=LineW, label=LegText[i])

    Ptitle.leftjust(Paxes)
    Xaxis.config_axis(Paxes)
    Yaxis.config_axis(Paxes)

    if (LegLoc != 'none'):
        Paxes.legend(LegText, loc=LegLoc)


# ###################################
# 
# function [] = PlotFsFigBgraph(Paxes, Y, Bcolors, Pmarker, Ptitle, Xlabel, Blabels, Ylabel, Ylim, Fsize, ShowX, ShowY, LegText, LegLoc)
# 
#   axes(Paxes);
# 
#   % Want to show accummulation of factors using a plot like a bar graph.
#   % Each factor is shown as a bar with it's particular height (positive
#   % or negative) going from the running sum of the previous factors. Like
#   % a stacked bar graph but the pieces of the bars are fanned out beside
#   % each other so that a negative piece doesn't obstruct a neighboring
#   % positive piece.
# 
#   % The Y matrix is two rows, with the starting height in the first
#   % row and the ending height in the second row. Ie, the number of
#   % bars we want is equal to the number of columns in Y.
#   % 
#   % Proportion the x-axis so that the bar widths are 0.8 and
#   % distances between bar centers is 1.
#   Nbars = size(Y,2);
#   Xmin = 0;
#   Xmax = Nbars + 1;
#   Xlim = [ Xmin Xmax ];
#   Xticks = [ 1:Nbars ];
# 
#   % Use patch to draw each bar at its specified height
#   hold on;
#   for i = 1:Nbars
#     Bcolor = str2rgb(Bcolors{i});
#     X1 = i - 0.4; 
#     X2 = i + 0.4; 
#     Y1 = Y(1,i);
#     Y2 = Y(2,i);
# 
#     Xbar = [ X1 X1 X2 X2 ];
#     Ybar = [ Y1 Y2 Y2 Y1 ];
# 
#     % Draw a bar for the factor magnitude
#     patch(Xbar,Ybar,Bcolor);
# 
# %    % On just the factors F1, F2, F12, superimpose
# %    % an arrow to denote the sign of the factor.
# %    if ((i ~= 1) && (i ~= Nbars))
# %      Xarrow1 = [  i  i ];
# %      Yarrow1 = [ Y1 Y2 ];
# %      Ainc = 0.08;
# %      % lines for the arrowhead
# %      if (Y1 <= Y2)
# %        % positive value, place arrowhead on top
# %        Xarrow2 = [  i-Ainc  i ];
# %        Yarrow2 = [ Y2-Ainc Y2 ];
# %
# %        Xarrow3 = [  i+Ainc  i ];
# %        Yarrow3 = [ Y2-Ainc Y2 ];
# %      else
# %        % negative value, place arrowhead on bottom
# %        Xarrow2 = [  i-Ainc  i ];
# %        Yarrow2 = [ Y2+Ainc Y2 ];
# %
# %       Xarrow3 = [  i+Ainc  i ];
# %       Yarrow3 = [ Y2+Ainc Y2 ];
# %     end
# % 
# %     line(Xarrow1, Yarrow1, 'Color', 'k');
# %     line(Xarrow2, Yarrow2, 'Color', 'k');
# %     line(Xarrow3, Yarrow3, 'Color', 'k');
# %    end
#   end
# 
#   set(Paxes, 'Xtick', Xticks);
#   set(Paxes, 'XTickLabel', Blabels);
#  
#   set(Paxes, 'FontSize', Fsize);
#   set(Paxes, 'LineWidth', 2);
#   set(Paxes, 'TickLength', [ 0.025 0.025 ]);
# 
#   xlabel(Xlabel);
#   xlim(Xlim);
# 
#   ylabel(Ylabel);
#   ylim(Ylim);
# 
#   if (isempty(Pmarker))
#     title(Ptitle);
#   else
#     T = title(sprintf('(%s) %s', Pmarker, Ptitle));
#     LeftJustTitle(Paxes, T);
#   end
# 
#   if (~strcmp(LegText, 'none'))
#     legend(LegText, 'Location', LegLoc);
#   end
# end
