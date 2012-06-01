function [ ] = PlotHist2d( Hist, Times, Bins, Clevs, Ptitle, OutFile )
%PlotHist2d create 2D contour plot of the 2D histograms
%   This routine will create 2D contour plots of the histograms created by GenHist2d.m.
%
%   The argument Times are the time step numbers which become the x values for the
%   2D plot.
%
%   The argument Bins defines the bin ranges for the histogram which become
%   The y values of the 2D plot.
%
%   The argument Ptitle is a title string for the plot.
%
%   The argument OutFile is the image file for the plot.
%

FontSize = 18;

% Show hours on the Time axis
%   Tick marks are too close when showing every hour, so show every 5 hours
%   Time increment is 5 min. so 5 hours is 60 steps
XtickLocs = find(mod(Times,60) == 0);
XtickVals = (XtickLocs - 1) / 12;


Fig = figure;
contourf(Times, Bins, Hist, Clevs), shading flat;
set(gca, 'XTick', XtickLocs);
set(gca, 'XTickLabel', XtickVals);
set(gca, 'FontSize', FontSize);
title(Ptitle);
xlabel('Sim Time (hrs)');
ylabel('Precip Rate (mm/hr)');
colorbar('FontSize', FontSize);

saveas(Fig, OutFile);

close(Fig);
