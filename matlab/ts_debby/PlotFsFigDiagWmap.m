function [] = PlotFsFigDiagWmap(Paxes, X, Y, Z, LatBounds, LonBounds, Clevs, Cmap, ClabInc, Pmarker, Ptitle, Fsize)

  axes(Paxes);

  LineW = 2;

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller', 'FontSize', 8);

  % Color levels and limits
  Cmin = Clevs(1);
  Cmax = Clevs(end);
  Clim = [ Cmin Cmax ];

  % replace small values with nans so that they don't appear
  % as dark regions in the contourfm plots
  Z(Z < Cmin) = nan;

  contourfm(Y, X, Z', Clevs, 'LineStyle', 'none'); 
  caxis(Paxes, Clim);

  % Put in a colorbar. The colorbar labeling has too many labels and comes out
  % cluttered. Just read in the Ytick and YtickLabel values that contourcmap
  % created, select every nth tick and label, and then reset the colorbar.
  Cbar = contourcmap(Cmap, 'ColorBar', 'on', 'Location', 'vertical');

  Cticks = get(Cbar, 'Ytick');
  CtickLabels = cellstr(get(Cbar, 'YtickLabel'))';

  Cticks = Cticks(1:ClabInc:end);
  CtickLabels = { CtickLabels{1:ClabInc:end} };

  set(Cbar, 'Ytick', Cticks);
  set(Cbar, 'YtickLabel', CtickLabels);

  % Draw the coast after the contourfm call so that coast lines will appear on top
  % of the contour plot
  Coast = load('coast');
  plotm(Coast.lat, Coast.long, 'Color', 'k', 'LineWidth', LineW);

  % Mark Africa
  textm(14, -12, 'Africa', 'FontSize', 10, 'Color', 'k', 'Rotation', 90);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end

