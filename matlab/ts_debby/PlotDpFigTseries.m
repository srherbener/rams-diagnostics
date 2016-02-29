function [] = PlotDpFigTseries(Paxes, X, Y, Pmarker, Ptitle, Ylabel, Fsize, ShowX, ShowY, LegText, LegLoc)

  axes(Paxes);

  plot(X, Y, 'LineWidth', 2);

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', 2);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);

  set(Paxes, 'XTick', [ 6 18 30 42 54 ]);
  if (ShowX > 0)
    set(Paxes, 'XTickLabel', { ' 12Z\newline22Aug' '  0Z\newline23Aug' ' 12Z\newline23Aug' '  0Z\newline24Aug' ' 12Z\newline24Aug' });
  else
    set(Paxes, 'XTickLabel', {});
  end

  ylabel(Ylabel);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

  if (~strcmp(LegLoc, 'none'))
    legend(LegText, 'Location', LegLoc);
  end

end
