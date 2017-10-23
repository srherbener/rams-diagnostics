function [ Fig Paxes ] = PlotDpFigMapContour(MapLoc, LAT, LON, Z, StormLoc, Pmarker, Ptitle, Cmin, Cinc, Cmax, CmapName, Fsize)

  Fig = figure;

  % contour levels
  Clevs = [ Cmin:Cinc:Cmax ];

  LatBounds = [ MapLoc(1) MapLoc(2) ];
  LonBounds = [ MapLoc(3) MapLoc(4) ];

  StormLat = [ StormLoc(1) ];
  StormLon = [ StormLoc(2) ];

  LineW = 3;
  Aratio = 1.4;

  % Location on page
  PlotLoc = [ 0.2 0.2 0.6 0.7 ];
  CbarLoc = [ 0.2 0.1 0.6 0.1 ];

  % Use cylindrical projection - Miller
  % Change aspect ratio of plotbox and data to make map longer in the x axes direction
  Paxes = worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');
  set(Paxes, 'Position', PlotLoc);
  set(Paxes, 'FontSize', Fsize);
  setm(Paxes, 'FontSize', Fsize);
  pbaspect( [ Aratio 1 1 ]);
  daspect( [ 1 Aratio 1 ]);

  % Plot on linear scale
  contourfm(LAT, LON, Z, Clevs, 'LineStyle', 'none'); 
  caxis([ Cmin Cmax ]);

  % Mark TS Debby location
  plotm(StormLat, StormLon, 'Color', 'r', 'LineStyle', 'none', 'Marker', 'x', 'LineWidth', LineW, 'MarkerSize', Fsize);

  % Colorbar, this sets the colormap too.
  Cbar = contourcmap('hot', 'ColorBar', 'on', 'Location', 'horizontal');
  set(Cbar, 'Position', CbarLoc);
  set(Cbar, 'FontSize', Fsize);
  xlabel(Cbar, 'M_D_S_F_C (kg km^-^2)');

%   % Reduce the amount of tick marks so that the labels won't overlap
%   Cticks = [ Cmin:Cinc:Cmax ];
%   CtickLabels = { '0.0' '0.6' '1.2' '1.8' '2.4' '3.0' };
%   set(Cbar, 'XTick', Cticks);
%   set(Cbar, 'XTickLabel', CtickLabels);

  if (strcmp(CmapName, 'hotrev'))
    % Reverse the 'hot' colormap for showing dust
    % Make sure this code comes after the contourcmap call since contourcmap
    % also sets the colormap.
    Cmap = colormap('hot');
    Cmap = Cmap(end:-1:1,:); % reverse the entries
  else
    Cmap = colormap(CmapName);
  end
  colormap(Paxes, Cmap);
  colormap(Cbar, Cmap);

  % title - left justify if a panel marker was specified
  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end
