function [] = PlotFsFigProfile(Paxes, X, Z, Pmarker, Ptitle, Xlabel, Xlim, Ylabel, Ylim, Fsize, LegText, Colors)

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

  legend(Plines, LegText,'Location', 'NorthEast', 'FontSize', LegendFsize);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
