function [ ] = PlotHist2d( Hist, Times, Bins, Ptitle, OutFile )
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



Fig = figure;
contourf(Times, Bins, Hist), shading flat;
title(Ptitle);
xlabel('Simulation Time');
ylabel('Precipitation Rate (mm/hr)');
colorbar;

saveas(Fig, OutFile);

close(Fig);
