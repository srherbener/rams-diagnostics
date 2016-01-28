function [] = PlotDpFigDustHov(Paxes, X, Y, Z, Pmarker, Ptitle, Fsize, ShowX, ShowY, Ylim)

  axes(Paxes);

  contourf(X, Y, Z, 20, 'LineStyle', 'none');

  % place colorbar beneath contour plot, but don't allow matlab to put in large gaps
  GapSize = 0.05;
  PaxesLoc = get(Paxes, 'Position'); % open up a small gap below plot
  CbarLoc = PaxesLoc;
  PaxesLoc(2) = PaxesLoc(2) + GapSize;
  PaxesLoc(4) = PaxesLoc(4) - GapSize;
  set(Paxes, 'Position', PaxesLoc);
  
  CbarLoc(4) = GapSize * 0.8;  % Force colorbar into the above gap
  Cbar = colorbar('Location', 'SouthOutside', 'Position', CbarLoc);
  caxis([ -3 0 ]);
  ylim(Ylim);

  set(Cbar, 'Ticks', [ -2 -1 0 ]);
  set(Cbar, 'TickLabels', { '10^-^2' '10^-^1' '1' });

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', 2);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);

  set(Paxes, 'XTick', [ 6 30 54 ]);
  if (ShowX > 0)
    set(Paxes, 'XTickLabel', { ' 12Z\newline22Aug' ' 12Z\newline23Aug' ' 12Z\newline24Aug' });
  else
    set(Paxes, 'XTickLabel', {});
  end

  set(Paxes, 'YTick', [ 0 5 10 15 ]);
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

end
