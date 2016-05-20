function [] = PlotFsFigXsection(Paxes, X, Y, Z, Pmarker, Ptitle, Xlabel, Xlim, Ylabel, Ylim, Cmap, Clim, Clevs, Fsize, ShowX, ShowY)

  axes(Paxes);

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', 2);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);

  contourf(X, Y, Z, Clevs, 'LineStyle', 'none');
  colorbar('EastOutside');
  caxis(Clim);

  if (strcmp(Cmap, 'redblue'))
    % need to call a function for redblue
    colormap(Paxes, redblue);
  else
    % expect a built in name
    colormap(Paxes, Cmap);
  end

  if (ShowX == 0)
    set(Paxes, 'XTickLabel', {});
  else
    xlabel(Xlabel);
  end
  xlim(Xlim);

  if (ShowY == 0)
    set(Paxes, 'YTickLabel', {});
  else
    ylabel(Ylabel);
  end
  ylim(Ylim);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end
end
