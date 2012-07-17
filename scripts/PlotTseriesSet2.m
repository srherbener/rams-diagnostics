function [ ] = PlotTseriesSet2( Times, Ts, Ptitle, Ylim, Ylabel, Lcolors, LegText, LegLoc, OutFile )
%PlotTseriesSet2 Plot a set of time series on the same panel
%   This function will take data contained in Ts and plot them on a single
%   panel.
%
%   Ccn and Sst define conditions (CCN concentrations and Sea Surface
%   Temperatures) associated with the histograms. These values will be used
%   for the plot legend.
%
%   Ptitle is a string (or array of strings) holding the title for the
%   plot.
%
%   OutFile is the path to the file that contains the image of the plot.
%

Fig = figure;

Lwidth = 2;
FontSz = 20;

Nts = size(Ts,1);

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first hist outside the loop, issue a "hold
% on" and then plot the remainder hists.

plot(Times, Ts(1,:), 'Color', char(Lcolors(1)), 'LineWidth', Lwidth);
set(gca, 'FontSize', FontSz);
ylim(Ylim);

hold on;

for i = 2:Nts % each column is a separate time series
    plot(Times, Ts(i,:), 'Color', char(Lcolors(i)),'LineWidth', Lwidth);
end

legend(LegText, 'Location', char(LegLoc));
legend boxoff;
title(Ptitle);
xlabel('Simulation Time (hr)');
ylabel(Ylabel);

saveas(Fig, OutFile);

hold off;
close(Fig);

end

