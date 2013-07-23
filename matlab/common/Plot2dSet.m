function [ ] = Plot2dSet( X, Y, Ptitle, Pmarkers, Xlabel, Ylabel, Lcolors, Lstyles, Gscales, LegText, LegLoc, AxisProps, AddMeas, OutFile )
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

Lwidth = 4;
LegFsize = 25;
Nprops = length(AxisProps);
Nplots = size(Y,1);

PanelTitle = ~isempty(Pmarkers);
if (PanelTitle)
    Ptitle = sprintf('(%s) %s', Pmarkers{1}, Ptitle);
end

% Need to establish an axis style before calling "hold on" (for the
% subsequent plots). Plot the first hist outside the loop, issue a "hold
% on" and then plot the remainder hists.

if (strcmp(Lcolors{1}, 'k'))
  Lcolor = [ 1 1 1 ] * Gscales(1);
else
  Lcolor = Lcolors{1};
end
plot(X(1,:), Y(1,:), 'Color', Lcolor, 'LineStyle', Lstyles{1}, 'LineWidth', Lwidth);
for i = 1:Nprops
  set(gca, AxisProps(i).Name, AxisProps(i).Val);
end
set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [ 0.025 0.025 ]);

hold on;

for i = 2:Nplots % each row is a separate curve for plotting
    if (strcmp(Lcolors{i}, 'k'))
      Lcolor = [ 1 1 1 ] * Gscales(i);
    else
      Lcolor = Lcolors{i};
    end
    plot(X(i,:), Y(i,:), 'Color', Lcolor, 'LineStyle', Lstyles{i}, 'LineWidth', Lwidth);
end

legend(LegText, 'Location', char(LegLoc), 'FontSize', LegFsize);
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

% Add in the temporal phases if requested
if (strcmp(AddMeas, 'Tphases'))
    Ylims = ylim;
    Yrange = Ylims(2) - Ylims(1);
    Yinc = 0.02 * Yrange;
    % put text opposite the place where the legend is located
    if (regexp(char(LegLoc), '[nN]orth'))
        % text goes on bottom
        Ty = Ylims(1) + 0.05 * Yrange;
    else
        % text goes on top
        Ty = Ylims(1) + 0.95 * Yrange;
    end
    
    %DrawTmark( 40,  60, 5, Ty, Yinc, 'RI');
    %DrawTmark( 80, 100, 5, Ty, Yinc, 'TR');
    DrawTmark(120, 140, 5, Ty, Yinc, 'SS');
end

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

%%%%%%%%%%%%%%%%%%%%
function [] = DrawTmark(X1, X2, Xinc, Y, Yinc, Text)

X  = (X1 + X2) / 2;
Y1 = Y - Yinc;
Y2 = Y + Yinc;

Bcolor = [ 1 1 1 ] * 0.25;

% Put text in middle
text(X, Y, Text, 'Color', Bcolor, 'FontSize', 20, ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

% Draw brackets around the text
line([ X1 X1 ], [ Y1 Y2 ], 'Color', Bcolor, 'LineWidth', 2, 'LineStyle', '-');
line([ X2 X2 ], [ Y1 Y2 ], 'Color', Bcolor, 'LineWidth', 2, 'LineStyle', '-');

line([ X1 X1+Xinc ], [ Y  Y  ], 'Color', Bcolor, 'LineWidth', 2, 'LineStyle', '-');
line([ X2-Xinc X2 ], [ Y  Y  ], 'Color', Bcolor, 'LineWidth', 2, 'LineStyle', '-');

end

