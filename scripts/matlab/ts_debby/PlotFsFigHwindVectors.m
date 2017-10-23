function [] = PlotFsFigHwindVectors(Paxes, X, Y, U, V, Pmarker, Ptitle, Fsize, Vscale)

  axes(Paxes);

  LatBounds = [ 5 26 ];
  LonBounds = [ -42 -8 ];

  LineW = 2;

  set(Paxes, 'FontSize', Fsize);
  worldmap(LatBounds, LonBounds);
  setm(Paxes, 'MapProjection', 'miller');
  load coast;  % creates vars "lat" and "long" that contain coast lines
  plotm(lat, long, 'Color', 'k', 'LineWidth', LineW);

  % quiverm wants Lat first, Lon second
  %    X is Lon
  %    Y is lat
  %    U is delta along Lon
  %    V is delta along Lat
  [ YGRID XGRID ] = meshgrid(Y, X);
  quiverm(double(YGRID), double(XGRID), double(V), double(U), Vscale);

  if (isempty(Pmarker))
    title(Ptitle);
  else
    T = title(sprintf('(%s) %s', Pmarker, Ptitle));
    LeftJustTitle(Paxes, T);
  end
end
