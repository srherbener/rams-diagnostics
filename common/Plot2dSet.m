function [ ] = Plot2dSet( X, Y, Ptitle, Xlabel, Ylabel, Lspecs, LegText, LegLoc, AxisProps, OutFile )
%Plot2dSet Plot a set of 2D line plots on the same panel
%   This function will take data contained in X and Y and plot them on a
%   single panel.
%
%   Ptitle is a string (or array of strings) holding the title for the
%   plot.
%
%   AxisProps is a structure contain a list of axis property names and
%   associated values that are desired to be set.
%
%   OutFile is the path to the file that contains the image of the plot.
%

Fig = figure;

Lwidth = 2;
Nprops = length(AxisProps);
Nplots = size(Y,1);

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first hist outside the loop, issue a "hold
% on" and then plot the remainder hists.

plot(X(1,:), Y(1,:), Lspecs{1}, 'LineWidth', Lwidth);
for i = 1:Nprops
  set(gca, AxisProps(i).Name, AxisProps(i).Val);
end

hold on;

for i = 2:Nplots % each row is a separate curve for plotting
    plot(X(i,:), Y(i,:), Lspecs{i}, 'LineWidth', Lwidth);
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

