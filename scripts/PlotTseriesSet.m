function [ ] = PlotTseriesSet( Ts, Ccn, Ptitle, Ylabel, Lstyles, LegLoc, OutFile )
%PlotTseriesSet Plot a set of time series on the same panel
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
%   Lstyles is a list of line styles to use for the traces in the plot.
%
%   OutFile is the path to the file that contains the image of the plot.
%

Fig = figure;

Times = (0:size(Ts,1)-1);
Lwidth = 2;

% Show hours on the Time axis
%   Tick marks are too close when showing every hour, so show every 5 hours
%   Time increment is 5 min. so 5 hours is 60 steps
XtickLocs = find(mod(Times,60) == 0);
XtickVals = (XtickLocs - 1) / 12;

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first hist outside the loop, issue a "hold
% on" and then plot the remainder hists.

plot(Times, Ts(:,1), char(Lstyles(1)), 'LineWidth', Lwidth);
set(gca, 'FontSize', 20);
set(gca, 'XTick', XtickLocs);
set(gca, 'XTickLabel', XtickVals);

hold on;

Ltext(1) = { sprintf('CCN: %d/cc', Ccn(1)) };

for i = 2:size(Ts,2) % each column is a separate time series
    plot(Times, Ts(:,i), char(Lstyles(i)),'LineWidth', Lwidth);
    Ltext(i) = { sprintf('CCN: %d/cc', Ccn(i)) };
end

legend(Ltext, 'Location', LegLoc);
legend boxoff;
title(Ptitle);
xlabel('Simulation Time (hr)');
ylabel(Ylabel);

saveas(Fig, OutFile);

hold off;
close(Fig);

end

