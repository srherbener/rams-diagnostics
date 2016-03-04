function [] = PlotDpFigDustHov(Paxes, X, Y, Z, Pmarker, Ptitle, Fsize, ShowX, ShowY, ShowLev, CaxisRange)

  axes(Paxes);

  if (CaxisRange <= 1)
    Z = log10(Z);   % convert to log values
  end
  contourf(X, Y, Z, 20, 'LineStyle', 'none');

%  % place colorbar beneath contour plot, but don't allow matlab to put in large gaps
%  GapSize = 0.05;
%  PaxesLoc = get(Paxes, 'Position'); % open up a small gap below plot
%  CbarLoc = PaxesLoc;
%  PaxesLoc(2) = PaxesLoc(2) + GapSize;
%  PaxesLoc(4) = PaxesLoc(4) - GapSize;
%  set(Paxes, 'Position', PaxesLoc);
%
%  CbarLoc(4) = GapSize * 0.8;  % Force colorbar into the above gap
%  Cbar = colorbar('Location', 'SouthOutside', 'Position', CbarLoc);
%  % CaxisRange

  % Place colorbar on right side
  Cbar = colorbar;

  %   0 : log scale 1e-3 to  1
  %   1 : log scale 1e-1 to 100
  %   2 : linear scale 0 500
  if (CaxisRange == 0)
    caxis([ -3 0 ]);
    set(Cbar, 'Ticks', [ -2 -1 0 ]);
    set(Cbar, 'TickLabels', { '10^-^2' '10^-^1' '1' });
  elseif (CaxisRange == 1)
    caxis([ -1 2 ]);
    set(Cbar, 'Ticks', [ 0 1 2 ]);
    set(Cbar, 'TickLabels', { '1' '10' '10^2' });
  else
    caxis([ 0 500 ]);
    set(Cbar, 'Ticks', [ 200 300 400 500 ]);
  end

  % ShowLev
  %   < 0  : 0  7   km
  %  == 0  : 0 16.5 km
  %   > 0  : 7 16.5 km
  if (ShowLev < 0)
    ylim([ 0 7 ]);
  elseif (ShowLev == 0)
    ylim([ 0 16.5 ]);
  else
    ylim([ 7 6.5 ]);
  end

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', 2);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);

  set(Paxes, 'XTick', [ 6 18 30 42 54 ]);
  if (ShowX > 0)
    set(Paxes, 'XTickLabel', { ' 12Z\newline22Aug' '  0Z\newline23Aug' ' 12Z\newline23Aug' '  0Z\newline24Aug' ' 12Z\newline24Aug' });
  else
    set(Paxes, 'XTickLabel', {});
  end

  if (ShowY > 0)
    ylabel('Height (km)');
  else
    set(Paxes, 'YTickLabel', {});
  end

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

  % Mark the separation level (7 km elevation)
  line([ 0 60 ], [ 7 7 ], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);

end
