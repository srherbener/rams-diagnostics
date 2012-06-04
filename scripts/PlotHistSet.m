function [ ] = PlotHistSet( Hists, Bins, Ccn, Sst, Npts, Lstyles, OutFile )
%PlotHistSet Plot a set of histograms on the same panel
%   This function will take data contained in Hist (y values) and Bins(x
%   values) and plot them on a single panel line plot. The Columns of Hist
%   and Bins are the separate histogram functions.
%
%   Ccn and Sst define conditions (CCN concentrations and Sea Surface
%   Temperatures) associated with the histograms. These values will be used
%   for the plot legend.
%
%   Npts contains the total number of points in the data set used to create
%   each histogram. This will be used to convert the histogram counts into
%   fractional area values.
%
%   OutFile is the path to the file that contains the image of the plot.
%
%   Lstyles is a list of line styles to use for the traces in the plot.

Fig = figure;

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first hist outside the loop, issue a "hold
% on" and then plot the remainder hists.

EndBin = 29; % just take the lower rain rates: 0 - 2.9 mm/hr

semilogy(Bins(1:EndBin,1), Hists(1:EndBin,1)/Npts(1), char(Lstyles(1)), 'LineWidth', 2);
set (gca, 'FontSize', 20);
hold on;
Ltext(1) = { sprintf('CCN: %d/cc', Ccn(1)) };

for i = 2:size(Hists,2) % each column is a separate histogram
    semilogy(Bins(1:EndBin,i), Hists(1:EndBin,i)/Npts(i), char(Lstyles(i)),'LineWidth', 2);
    Ltext(i) = { sprintf('CCN: %d/cc', Ccn(i)) };
end

legend(Ltext);
Ptitle = sprintf('Rainfall Rate PDF, SST: %d K',Sst(1));
title(Ptitle);
xlabel('Rainfall Rate (mm/hr)');
ylabel('Fraction of domain');

saveas(Fig, OutFile);

hold off;
close(Fig);

end

