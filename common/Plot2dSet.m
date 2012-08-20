function [ ] = Plot2dSet( X, Y, Ptitle, Xlabel, Ylabel, Lcolors, LegText, LegLoc, OutFile )
%Plot2dSet Plot a set of 2D line plots on the same panel
%   This function will take data contained in X and Y and plot them on a
%   single panel.
%
%   Ptitle is a string (or array of strings) holding the title for the
%   plot.
%
%   OutFile is the path to the file that contains the image of the plot.
%

Fig = figure;

Lwidth = 2;
FontSz = 20;

Nplots = size(Y,1);

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first hist outside the loop, issue a "hold
% on" and then plot the remainder hists.

plot(X(1,:), Y(1,:), 'Color', char(Lcolors(1)), 'LineWidth', Lwidth);
set(gca, 'FontSize', FontSz);

hold on;

for i = 2:Nplots % each row is a separate curve for plotting
    plot(X(i,:), Y(i,:), 'Color', char(Lcolors(i)),'LineWidth', Lwidth);
end

legend(LegText, 'Location', char(LegLoc));
legend boxoff;
title(Ptitle);
xlabel(Xlabel);
ylabel(Ylabel);

saveas(Fig, OutFile);

hold off;
close(Fig);

end

