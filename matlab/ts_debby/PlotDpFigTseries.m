function [] = PlotDpFigTseries(Paxes, X, Y, Lcolors, Pmarker, Ptitle, Ylabel, Yscale, Ylim, Fsize, ShowX, ShowY, LegText, LegLoc)

  axes(Paxes);

  Ysize = size(Y);
  Nlines = Ysize(2);

  set(Paxes, 'Yscale', Yscale);
  hold on;
  for iline = 1:Nlines
    line(X, Y(:,iline), 'LineWidth', 2, 'Color', str2rgb(Lcolors{iline}));
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

  if (strcmp(Yscale,'log'))
    set(Paxes, 'Ytick', [ 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 ]);
  end

  ylabel(Ylabel);
  ylim(Ylim);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

  if (~strcmp(LegLoc, 'none'))
    Leg = legend(LegText, 'Location', LegLoc);
  end

  hold off;
end
