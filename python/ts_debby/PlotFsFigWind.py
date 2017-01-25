#!/usr/bin/env python3

import sys
sys.path.append("/Users/steve/etc/python/common")

import matplotlib.pyplot as plt
import numpy as np
import PlotUtils as plu
import h5py

CaseList = [ [ 'TSD_SAL_DUST',      'SD',   'black' ],
             [ 'TSD_SAL_NODUST',    'SND',  'blue'  ],
             [ 'TSD_NONSAL_DUST',   'NSD',  'red'   ],
             [ 'TSD_NONSAL_NODUST', 'NSND', 'green' ] ]
Ncases = len(CaseList)

LegText = [ ]
Colors = [ ]
for i in range(0, Ncases):
    Case  = CaseList[i][0]
    Label = CaseList[i][1]
    Color = CaseList[i][2]

    InFname = "DIAGS/storm_meas_{0:s}.h5".format(Case)
    InVname = "/max_wind_t"
    print("Reading {0:s} ({1:s})".format(InFname, InVname))

    InFile = h5py.File(InFname, mode='r')

    # grab the time coordinate values, and initialize the wind array
    if (i == 0):
        T = InFile['/t_coords'][:] / 3600 - 42  # sim time in hours starting with zero
        Nt = len(T)
        WindTs = np.zeros( [ Nt, Ncases ], dtype=np.float_)

    WindTs[ :, i ] = InFile[InVname][:]
    LegText.append(Label)
    Colors.append(Color)

print("")

# Create 3 panel plot
Pdir = 'Plots.py'

Fig = plt.figure()

# Vt max time series, top half
Paxes = Fig.add_subplot(2,1,1)

Xaxis = plu.AxisConfig('x', [ T[0], T[-1] ], r'$Sim Time (h)$')
Yaxis = plu.AxisConfig('y', [ 0, 25 ], r'$Speed (ms^{-1})$')
Ptitle = plu.TitleConfig('a', r'$V_t (max)$')

Fsize = 20

LegLoc = 'lower right'
plu.PlotLine(Paxes, T, WindTs, Ptitle, Xaxis, Yaxis, Fsize, LegText, LegLoc, Colors)


OutFile = "{0:s}/FsFigWindMax.png".format(Pdir)
print("Writing: {0:s}".format(OutFile))
Fig.savefig(OutFile)

# 
#   Ncases = length(CaseList);
#   
#   Fsize = 13;
# 
#   WmaxMeas = 'max';  % 'avg' -> use avg Vt measurement
#                      % 'max' -> use max Vt measurement
# 
#   % Read in the max wind time series
#   for icase = 1:Ncases
#     Case  = CaseList{icase}{1};
#     Label = CaseList{icase}{2};
#     Color = CaseList{icase}{3};
# 
#     InFname = sprintf('DIAGS/storm_meas_%s.h5', Case);
#     InVname = sprintf('/%s_wind_t', WmaxMeas);
#     fprintf('Reading %s (%s)\n', InFname, InVname);
#     if (icase == 1)
#       T = squeeze(h5read(InFname, '/t_coords'))./3600 - 42; % convert to sim time in hours starting with zero
#       Nt = length(T);
#       WIND_TS = zeros([Nt Ncases]);
#     end
#     WIND_TS(:,icase)   = squeeze(h5read(InFname, InVname));
#     WmaxLegText{icase} = Label;
#     WmaxColors{icase}  = Color;
#   end
#   WmaxLegLoc  = 'NorthEastOutside';
#   fprintf('\n');
# 
#   % Read in factor separation data
#   InFname = 'DIAGS/fs_factors.h5';
# 
#   InVname = sprintf('/s_%s_wind_t_bar_factors', WmaxMeas);
#   fprintf('Reading %s (%s)\n', InFname, InVname);
#   SAL_FS_WIND = squeeze(h5read(InFname, InVname));
# 
#   InVname = sprintf('/ps_%s_wind_t_bar_factors', WmaxMeas);
#   fprintf('Reading %s (%s)\n', InFname, InVname);
#   PRESAL_FS_WIND = squeeze(h5read(InFname, InVname));
# 
#   % Plot: 2 panels (2x1)
#   Fig = figure;
# 
#   % pull up the bottom of the top panel and push down the top
#   % of the bottom panel in order to make room for labeling
#   PlocInc = 0.01;
#   switch(WmaxMeas)
#     case 'avg'
#       OutFile = sprintf('%s/FsFigWindAvg.jpg', Pdir);
#       WmaxYlim = [ 0 15 ];
#     case 'max'
#       OutFile = sprintf('%s/FsFigWindMax.jpg', Pdir);
#       WmaxYlim = [ 0 25 ];
#   end
#   Ptitle = sprintf('V_t (%s)', WmaxMeas);
# 
#   Paxes = subplot(2, 2, [ 1 2 ]);
#   Ploc = get(Paxes, 'Position');
#   Ploc(2) = Ploc(2) + PlocInc;
#   Ploc(3) = Ploc(3) * 0.97;    % get the ends of the plots to line up
#   Ploc(4) = Ploc(4) - PlocInc;
#   set(Paxes, 'Position', Ploc);
#   PlotFsFigTseries(Paxes, T, WIND_TS, WmaxColors, 'a', Ptitle, 'Speed (ms^-^1)', 'linear', WmaxYlim, Fsize, 1, 1, WmaxLegText, WmaxLegLoc);
# 
#   % Mark the PRESAL (10-30 h) and SAL (40-60 h) time periods
#   line([ 10 30 ], [ 5 5 ], 'Color', 'k', 'LineWidth', 2);
#   line([ 10 10 ], [ 4 6 ], 'Color', 'k', 'LineWidth', 2);
#   line([ 30 30 ], [ 4 6 ], 'Color', 'k', 'LineWidth', 2);
#   text(15, 6.5, 'PSAP', 'FontSize', Fsize);
# 
#   line([ 40 60 ], [ 5 5 ], 'Color', 'k', 'LineWidth', 2);
#   line([ 40 40 ], [ 4 6 ], 'Color', 'k', 'LineWidth', 2);
#   line([ 60 60 ], [ 4 6 ], 'Color', 'k', 'LineWidth', 2);
#   text(48, 6.5, 'SAP', 'FontSize', Fsize);
# 
#   % Bar graphs
#   BarColors = { 'dodgerblue' 'cyan' 'orange' 'yellow' 'red' };
#   BarLabels = { 'NSND' 'F1' 'F2' 'F12' 'SD' };
#   BarLegText = 'none';
#   BarLegLoc = '';
#   switch(WmaxMeas)
#     case 'avg'
#       Bylim = [ 9 11 ];
#     case 'max'
#       Bylim = [ 13 18 ];
#   end
# 
#   Paxes = subplot(2, 2, 3);
#   Ploc = get(Paxes, 'Position');
#   Ploc(4) = Ploc(4) - PlocInc;
#   set(Paxes, 'Position', Ploc);
#   PlotFsFigBgraph(Paxes, PRESAL_FS_WIND, BarColors, 'b', 'PSAP', 'Factor', BarLabels, 'Magnitude', Bylim, Fsize, 1, 1, BarLegText, BarLegLoc );
# 
#   Paxes = subplot(2, 2, 4);
#   Ploc = get(Paxes, 'Position');
#   Ploc(4) = Ploc(4) - PlocInc;
#   set(Paxes, 'Position', Ploc);
#   PlotFsFigBgraph(Paxes, SAL_FS_WIND, BarColors, 'c', 'SAP', 'Factor', BarLabels, 'Magnitude', Bylim, Fsize, 1, 1, BarLegText, BarLegLoc );
# 
#   fprintf('Writing: %s\n', OutFile);
#   saveas(Fig, OutFile);
#   close(Fig);
# end
