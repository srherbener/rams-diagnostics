function [ ] = PlotContour( X, Y, Z, Ptitle, Pmarker, Xlabel, Ylabel, Fill, Cbar, Cmap, Crange, Clevs, AxisProps, OutFile )
%PlotContour Creae a contour plot
%   This function will plot data contained in Z as contours on a grid defined
%   by X and Y.
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

if (strcmp(Cmap, 'LightGray'))
  % Gray scale, but go from [ g g g ] to [ 1 1 1 ] instead
  % of [ 0 0 0 ] to [ 1 1 1 ] where g > 0.
  % This will get rid of dark regions.
  %
  % x = 1 to 64, corresponds to y = g to 1.0
  % produce a linear scale
  %   y = mx + b
  %       m = (y2 - y1)/(x2 - x1)
  %       b = y1 - m*x1

  x1 = 1;
  x2 = 64;
  y1 = 0.3; % <-- g  (0.3 seems to get good contrast without making the overall picture too dark)
  y2 = 1.0;

  m = (y2 - y1) / (x2 - x1);
  b = y1 - (m * x1);
  
  x = x1:x2;
  y = m .* x + b;

  ColorMap = repmat(y', [ 1 3 ]);
else
  ColorMap = Cmap;
end

Nprops = length(AxisProps);

PanelTitle = ~isempty(Pmarker);
if (PanelTitle)
    Ptitle = sprintf('(%s) %s', Pmarker, Ptitle);
end

% create the plot
if (Fill == 1)
  % filled contours
  contourf(X, Y, Z, Clevs, 'LineColor', 'none');
end

if (~isempty(Crange))
  caxis(Crange);
end

if (Cbar == 1)
  colorbar;
end

% Set up the axis
for i = 1:Nprops
  set(gca, AxisProps(i).Name, AxisProps(i).Val);
end

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

% If panel style of title, Fix up the positioning
if (PanelTitle)
  Ppos = get(gca, 'Position'); % position of plot area
  Ppos(1) = Ppos(1) * 1.00;
  Ppos(2) = Ppos(2) * 1.00;
  Ppos(3) = Ppos(3) * 0.90;
  Ppos(4) = Ppos(4) * 0.90;
  set(gca, 'Position', Ppos);
end

saveas(Fig, OutFile);
close(Fig);

end
