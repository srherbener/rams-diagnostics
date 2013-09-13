function [ ] = PlotBarSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Bcolors, LegText, LegLoc, AxisProps, OutFile )
%PlotBarSet Plot a set of bar plots on the same panel
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

Nprops = length(AxisProps);
Nplots = size(Y,1);

PanelTitle = ~isempty(Pmarkers);
if (PanelTitle)
    Ptitle = sprintf('(%s) %s', Pmarkers{1}, Ptitle);
end

if (~isempty(Bcolors)
  bar(X, Y, 'FaceColor', Bcolors)
end

for i = 1:Nprops
  set(gca, AxisProps(i).Name, AxisProps(i).Val);
end

legend(LegText, 'Location', char(LegLoc), 'FontSize', 25);
legend boxoff;

if (~strcmp(Ptitle, ' '))
  if (PanelTitle)
      % The title is in a box that adjusts to the amount of characters in
      % the title. Ie, it doesn't do any good to do Left/Center/Right
      % alignment. But, the entire box can be moved to the left side of the
      % plot.
      T = title(Ptitle);
      set(T, 'Units', 'Normalized');
      set(T, 'HorizontalAlignment', 'Left');
      Tpos = get(T, 'Position');
      Tpos(1) = 0; % line up with left edge of plot area
      set(T, 'Position', Tpos);
  else
      title(Ptitle);
  end
end
xlabel(Xlabel);
ylabel(Ylabel);

% Fix up the positioning
Ppos = get(gca, 'Position'); % position of plot area
Ppos(1) = Ppos(1) * 1.00;
Ppos(2) = Ppos(2) * 1.00;
Ppos(3) = Ppos(3) * 0.90;
Ppos(4) = Ppos(4) * 0.90;
set(gca, 'Position', Ppos);

saveas(Fig, OutFile);

hold off;
close(Fig);

end
