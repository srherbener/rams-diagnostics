function [] = PlotFsFigProfile(Paxes, X, Z, Pmarker, Ptitle, Xlabel, Xlim, Ylabel, Ylim, Fsize, LegText, LegLoc, Colors)

  axes(Paxes);

  LegendFsize = 8;
  LineW = 2;

  [ Npts Nlines ] = size(X);

  set(Paxes, 'FontSize', Fsize);
  for i = 1:Nlines
    Xline = squeeze(X(:,i));
    Color = str2rgb(Colors{i});
    Plines(i) = line(Xline, Z, 'Color', Color, 'LineWidth', LineW);
  end

  xlim(Xlim);
  ylim(Ylim);

  xlabel(Xlabel);
  ylabel(Ylabel);

  if (~strcmp(LegLoc, 'none'))
    legend(Plines, LegText,'Location', LegLoc, 'FontSize', LegendFsize);
  end

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
