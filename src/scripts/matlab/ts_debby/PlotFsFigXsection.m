function [] = PlotFsFigXsection(Paxes, X, Y, Z, Pmarker, Ptitle, Xlabel, Xlim, Ylabel, Ylim, Cmap, Clim, Clevs, Fsize, ShowX, ShowY)

  axes(Paxes);

  set(Paxes, 'FontSize', Fsize);
  set(Paxes, 'LineWidth', 2);
  set(Paxes, 'TickLength', [ 0.025 0.025 ]);

  contourf(X, Y, Z, Clevs, 'LineStyle', 'none');
  colorbar('EastOutside', 'FontSize', Fsize);
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
    xlabel(Xlabel, 'FontSize', Fsize);
  end
  if (Xlim(1) > Xlim(2))
    set(Paxes, 'Xdir', 'reverse');
    Xlim = [ Xlim(2) Xlim(1) ];
  end
  xlim(Xlim);

  if (ShowY == 0)
    set(Paxes, 'YTickLabel', {});
  else
    ylabel(Ylabel, 'FontSize', Fsize);
  end
  if (Ylim(1) > Ylim(2))
    set(Paxes, 'Ydir', 'reverse');
    Ylim = [ Ylim(2) Ylim(1) ];
  end
  ylim(Ylim);

  if (isempty(Pmarker))
    title(Ptitle, 'FontSize', Fsize);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle), 'FontSize', Fsize);
    LeftJustTitle(Paxes, T);
  end
end
