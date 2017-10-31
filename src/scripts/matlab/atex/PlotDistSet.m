function [ ] = PlotDistSet( Hists, Bins, Ptitle, Xlabel, Ylabel, Ylim, Lcolors, Ltext, Lloc, OutFile )
%PlotDistSet Plot a set of histograms on the same panel
%   This function will take data contained in Hist (y values) and Bins(x
%   values) and plot them on a single panel line plot. The Columns of Hist
%   and Bins are the separate histogram functions.
%
%   Ccn and Sst define conditions (CCN concentrations and Sea Surface
%   Temperatures) associated with the histograms. These values will be used
%   for the plot legend.
%
%   OutFile is the path to the file that contains the image of the plot.
%

Fig = figure;

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first hist outside the loop, issue a "hold
% on" and then plot the remainder hists.

semilogy(Bins(:,1), Hists(:,1), 'Color', char(Lcolors(1)) , 'LineWidth', 2);
set (gca, 'FontSize', 20);
hold on;

for i = 2:size(Hists,2) % each column is a separate histogram
    semilogy(Bins(:,i), Hists(:,i), 'Color', char(Lcolors(i)),'LineWidth', 2);
end

legend(Ltext);
title(Ptitle);
xlabel(Xlabel);
ylabel(Ylabel);
ylim(Ylim);

saveas(Fig, OutFile);

hold off;
close(Fig);

end

