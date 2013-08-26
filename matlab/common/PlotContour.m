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

Nprops = length(AxisProps);

PanelTitle = ~isempty(Pmarker);
if (PanelTitle)
    Ptitle = sprintf('(%s) %s', Pmarker, Ptitle);
end

% create the plot
if (Fill == 1)
  % filled contours
  contourf(X, Y, Z, Clevs);
  shading flat;
else
  % line contours
  contour(X, Y, Z, Clevs);
end

if (~isempty(Cmap))
  colormap(Cmap);
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
