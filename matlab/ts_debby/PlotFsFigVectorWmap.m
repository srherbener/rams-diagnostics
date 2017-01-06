function [] = PlotFsFigVectorWmap(Paxes, Lat, Lon, LatComp, LonComp, LatBounds, LonBounds, Pmarker, Ptitle, Fsize, Qscale)

  axes(Paxes);

  LineW = 2;

  % grid for m_quiver
  [ Mlon Mlat ] = meshgrid(Lon, Lat);

  % draw the vectors on a rectangular map
  hold on;
  m_proj('miller', 'lat', LatBounds, 'lon', LonBounds);
  m_coast('color', 'k');
  m_grid('linestyle', 'none', 'box', 'fancy', 'tickdir', 'out');
  m_text(-14.7, 14, 'Africa', 'color', 'k', 'rotation', 90);
  m_quiver(Mlon, Mlat, LonComp', LatComp', Qscale); % Note transpose on *Comp vars
  hold off;

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(T);
  end

end
